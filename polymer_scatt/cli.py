from __future__ import annotations
import argparse
from collections import OrderedDict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import sys
import numpy as np
import scipy
from . import __version__
from .core import analyze_curve
from .io import iter_curves, file_hash
from .llm import interpret, build_prompt
from .report import summarize, group_name, safe_json, write_csv, write_html

def analyze(args):
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    results=[];representatives=OrderedDict();sources=[]
    for name in args.inputs:
        path=Path(name)
        sources.append({'filename':path.name,'bytes':path.stat().st_size,'sha256':file_hash(path)})
        for curve in iter_curves(path,q_unit=args.q_unit):
            result=analyze_curve(curve['q'],curve['intensity'],curve.get('sigma'),
                q_unit=curve['q_unit'],label=curve['label'],qmin=args.qmin,qmax=args.qmax)
            result['metadata']=curve['metadata']
            result=safe_json(result);results.append(result)
            group=group_name(result)
            if group not in representatives: representatives[group]=[(curve,result)]
            elif len(representatives[group])==1: representatives[group].append((curve,result))
            else: representatives[group][-1]=(curve,result)
            if len(results)%250==0: print(f'Analyzed {len(results)} curves',flush=True)
    if not results: raise ValueError('No supported curves found')
    summary=summarize(results)
    manifest={'software':'PolymerScattLab','version':__version__,
        'created_utc':datetime.now(timezone.utc).isoformat(),'inputs':sources,
        'configuration':{'q_unit_override':args.q_unit,'qmin':args.qmin,'qmax':args.qmax,
            'background_subtraction':'none added; use supplied reduced data',
            'third_column':'treated as independent 1-sigma I error; needs metadata confirmation',
            'frame_averaging':False,'sample_metadata':'not supplied'},
        'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'sources_modified':False,'external_llm_executed':False}
    interpretation=interpret(summary)
    (out/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False),encoding='utf-8')
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    (out/'llm_prompt.txt').write_text(json.dumps(build_prompt(summary),indent=2),encoding='utf-8')
    (out/'interpretation.json').write_text(json.dumps(interpretation,indent=2),encoding='utf-8')
    write_csv(results,out/'results.csv')
    reps=[p for pairs in representatives.values() for p in pairs]
    write_html(results,summary,reps,manifest,out/'report.html',interpretation)
    print(json.dumps(summary,indent=2));print(f'Report: {out.resolve()/"report.html"}')

def main(argv=None):
    p=argparse.ArgumentParser(description='Polymer SAXS/WAXS analysis with optional LLM interpretation')
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('analyze',help='Analyze DAT/CSV/XLSX/ZIP inputs and write CSV, JSON, HTML')
    a.add_argument('inputs',nargs='+');a.add_argument('--output','-o',default='analysis_output')
    a.add_argument('--q-unit',choices=['unknown','A^-1','nm^-1'],default='unknown',help='Only set from confirmed metadata; workbook headers take precedence')
    a.add_argument('--qmin',type=float);a.add_argument('--qmax',type=float)
    a.set_defaults(func=analyze)
    s=sub.add_parser('serve',help='Start local browser interface')
    s.add_argument('--port',type=int,default=8765);s.add_argument('--no-browser',action='store_true')
    s.set_defaults(func=lambda a:__import__('polymer_scatt.server',fromlist=['serve']).serve(port=a.port,open_browser=not a.no_browser))
    f=sub.add_parser('features',help='Explicit-range Debye hypothesis fitting or descriptive peak/d-spacing analysis')
    f.add_argument('inputs',nargs='+');f.add_argument('--q-unit',choices=['unknown','A^-1','nm^-1'],default='unknown')
    f.add_argument('--qmin',type=float,required=True);f.add_argument('--qmax',type=float,required=True)
    f.add_argument('--debye',action='store_true');f.add_argument('--fit-background',action='store_true')
    f.add_argument('--peaks',action='store_true');f.add_argument('--output',default='features.json')
    def features(a):
        from .core import fit_debye,find_peaks
        if not a.debye and not a.peaks: raise ValueError('Select --debye and/or --peaks')
        if not np.isfinite(a.qmin) or not np.isfinite(a.qmax) or a.qmin>=a.qmax: raise ValueError('Supply finite qmin < qmax')
        rows=[]
        for name in a.inputs:
            for c in iter_curves(name,q_unit=a.q_unit):
                r={'label':c['label'],'metadata':c['metadata'],'q_unit':c['q_unit'],
                   'user_selected_qmin':a.qmin,'user_selected_qmax':a.qmax}
                if a.debye: r['debye']=fit_debye(c['q'],c['intensity'],c['sigma'],c['q_unit'],a.qmin,a.qmax,a.fit_background)
                if a.peaks:
                    if c['q_unit']=='unknown': raise ValueError('Confirm q units before a physical d-spacing analysis')
                    sel=(c['q']>=a.qmin)&(c['q']<=a.qmax)
                    r['peaks']=find_peaks(c['q'][sel],c['intensity'][sel],c['q_unit'],
                        sigma=None if c['sigma'] is None else c['sigma'][sel])
                rows.append(safe_json(r))
        Path(a.output).write_text(json.dumps(rows,indent=2,allow_nan=False),encoding='utf-8')
        print(f'Wrote {len(rows)} exploratory records to {a.output}; no phase or crystallinity assignment.')
    f.set_defaults(func=features)
    e=sub.add_parser('explain',help='Optional LLM review of a summary.json; remote use must be explicit')
    e.add_argument('summary');e.add_argument('--provider',choices=['offline','ollama','openai-compatible'],default='offline')
    e.add_argument('--model');e.add_argument('--endpoint');e.add_argument('--question')
    e.add_argument('--api-key-env',default='POLYMERSCATT_API_KEY')
    e.add_argument('--output',default='llm_interpretation.json')
    def explain(a):
        summary=json.loads(Path(a.summary).read_text(encoding='utf-8'))
        r=interpret(summary,provider=a.provider,model=a.model,endpoint=a.endpoint,
            api_key=os.environ.get(a.api_key_env),question=a.question)
        Path(a.output).write_text(json.dumps(r,indent=2),encoding='utf-8');print(r['narrative'])
    e.set_defaults(func=explain)
    args=p.parse_args(argv)
    try: args.func(args)
    except (ValueError,OSError,RuntimeError) as exc:
        print(f'Error: {exc}',file=sys.stderr);return 1
    return 0

if __name__=='__main__': sys.exit(main())
