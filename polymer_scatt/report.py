"""Standalone reports; never substitute an LLM narrative for numerical results."""
from __future__ import annotations
import base64
from collections import Counter, defaultdict
import csv
import html
import io
import json
from pathlib import Path
import numpy as np

SOURCES = [
    ('SasView: Guinier model and scattering-density Rg', 'https://www.sasview.org/docs/user/models/guinier.html'),
    ('SasView: Debye Gaussian polymer coil model', 'https://www.sasview.org/docs/user/models/mono_gauss_coil.html'),
    ('BioXTAS RAW: Guinier validity and residual interpretation', 'https://bioxtas-raw.readthedocs.io/en/v2.1.4/saxs/saxs_guinier.html'),
    ('SasView: polydisperse Gaussian coils', 'https://www.sasview.org/docs/user/models/poly_gauss_coil.html'),
]

def safe_json(value):
    if isinstance(value, dict): return {str(k):safe_json(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [safe_json(v) for v in value]
    if isinstance(value,np.ndarray): return safe_json(value.tolist())
    if isinstance(value,(np.integer,)): return int(value)
    if isinstance(value,(float,np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value,(np.bool_,)): return bool(value)
    return value

def group_name(result):
    m=result.get('metadata',{})
    if m.get('sample') is not None: return f"Sample {m['sample']}"
    return str(result.get('label','Curve'))

def summarize(results):
    groups=defaultdict(list)
    for r in results: groups[group_name(r)].append(r)
    out=[]
    for name,rs in groups.items():
        reasons=Counter(reason for r in rs for reason in r.get('rejection_reasons',[]))
        out.append(dict(label=name,n_curves=len(rs),n_accepted=sum(r['status']=='accepted' for r in rs),
            n_rejected=sum(r['status']=='rejected' for r in rs),n_invalid=sum(r['status']=='invalid' for r in rs),
            q_unit=rs[0]['q_unit'],rejection_reasons=[f'{k} ({v} curves)' for k,v in reasons.most_common(5)]))
    return {'title':'Polymer SAXS/WAXS analysis','n_curves':len(results),
        'n_accepted':sum(r['status']=='accepted' for r in results),
        'n_rejected':sum(r['status']=='rejected' for r in results),
        'n_invalid':sum(r['status']=='invalid' for r in results),
        'groups':out,'sample_state':'unknown','chain_rg_status':'not established',
        'warnings':['Acceptance means an apparent Guinier fit passes engineering checks; molecular-chain assignment requires physical metadata.',
        'No frame averaging; frame numbers and T labels have no assumed time or temperature units.',
        'Confirm any unknown q units and the interpretation of supplied errors as independent 1-sigma intensity uncertainties.',
        'No absolute intensity calibration, concentration correction or additional background subtraction was applied.']}

def write_csv(results,path):
    fields=['label','group','source','member','frame','status','q_unit','rg_native','rg_nm','rg_stderr_native',
            'candidate_rg_native_NOT_VALIDATED','candidate_rg_nm_NOT_VALIDATED','qmin','qmax','qmin_rg','qmax_rg',
            'n_points','r2','reduced_chi2','endpoint_relative_spread','curvature_z','rejection_reasons','warnings']
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in results:
            c=r.get('fit') or r.get('candidate') or {};m=r.get('metadata',{})
            row={k:r.get(k) for k in ['label','status','q_unit','rg_native','rg_nm','rg_stderr_native']}
            row.update(group=group_name(r),source=m.get('source'),member=m.get('member'),frame=m.get('frame'))
            row.update({k:c.get(k) for k in ['qmin','qmax','qmin_rg','qmax_rg','n_points','r2','reduced_chi2','endpoint_relative_spread','curvature_z']})
            if r['status']!='accepted':
                row['candidate_rg_native_NOT_VALIDATED']=c.get('rg_native')
                row['candidate_rg_nm_NOT_VALIDATED']=c.get('rg_nm')
            row['rejection_reasons']='; '.join(r.get('rejection_reasons',[]))
            row['warnings']='; '.join(r.get('warnings',[]))
            # Avoid spreadsheet formula execution from uploaded labels/metadata.
            for k,v in row.items():
                if isinstance(v,str) and v.startswith(('=','+','-','@')): row[k]="'"+v
            w.writerow(row)

def diagnostic_png(curve,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    q=np.asarray(curve['q']);y=np.asarray(curve['intensity'])
    use=np.isfinite(q)&np.isfinite(y)&(q>0)&(y>0)
    q,y=q[use],y[use]
    fig,axs=plt.subplots(2,2,figsize=(11,7.4),layout='constrained')
    color='#087e8b';bad='#be5a24';unit=result['q_unit'];c=result.get('fit') or result.get('candidate')
    axs[0,0].loglog(q,y,color=color,lw=1)
    axs[0,0].set(xlabel=f'q ({unit})',ylabel='I(q), supplied units',title='Measured profile (positive intensities shown)')
    # Showing the true low-q region is essential; never hide it behind a preferred fit.
    limit=max(q[0]*4,1.05*c['qmax'] if c else 0) if len(q) else 0
    low=q<=min(q[-1],limit) if len(q) else np.zeros(0,dtype=bool)
    axs[0,1].plot(q[low]**2,np.log(y[low]),'.',ms=3,color=color)
    axs[0,1].set(xlabel=f'q² ({unit})²',ylabel='ln I(q)',title='Low-q Guinier diagnostic')
    axs[1,0].plot(q,q*q*y,color=color,lw=1)
    axs[1,0].set(xlabel=f'q ({unit})',ylabel='q² I(q), unnormalized',title='Kratky diagnostic')
    if c and c.get('slope') is not None and c.get('intercept') is not None:
        sel=(q>=c['qmin'])&(q<=c['qmax'])
        pred=c['intercept']+c['slope']*q[sel]**2
        axs[0,1].plot(q[sel]**2,pred,'--',color=bad,label='candidate' if result['status']!='accepted' else 'accepted apparent fit')
        axs[0,1].legend(fontsize=8)
        axs[1,1].plot(q[sel]**2,np.log(y[sel])-pred,'o-',ms=3,color=bad)
        axs[1,1].axhline(0,color='#777',lw=.7)
        axs[1,1].set(xlabel=f'q² ({unit})²',ylabel='ln I − fitted ln I',title='Fit residuals')
    else:
        axs[1,1].text(.5,.5,'No decreasing, usable Guinier candidate',ha='center',va='center',wrap=True,transform=axs[1,1].transAxes)
        axs[1,1].set_axis_off()
    for ax in axs.flat: ax.grid(alpha=.18)
    masked=result.get('diagnostics',{}).get('zero_zero_masked_bins',0)
    nonpositive=result.get('diagnostics',{}).get('nonpositive_intensity_bins',0)
    fig.suptitle(f"{result['label']}  |  {result['status'].upper()}  |  masked: {masked}; other nonpositive: {nonpositive}",fontsize=12,fontweight='bold')
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=125);plt.close(fig)
    return buf.getvalue()

def fmt(x,digits=4):
    return '—' if x is None else f'{x:.{digits}g}' if isinstance(x,(int,float)) else str(x)

def write_html(results,summary,reps,manifest,path,interpretation=None,notes=''):
    esc=lambda s:html.escape(str(s))
    counts=summary
    rows=''.join('<tr>'+''.join(f'<td>{esc(g[k])}</td>' for k in ['label','n_curves','n_accepted','n_rejected','n_invalid','q_unit'])+'</tr>' for g in summary['groups'])
    rep_sections=[]
    for curve,r in reps:
        c=r.get('fit') or r.get('candidate') or {}
        img=base64.b64encode(diagnostic_png(curve,r)).decode()
        nm=c.get('rg_nm');unit=r.get('length_unit','inverse native q')
        rgtext=f"{fmt(c.get('rg_native'))} {unit}"+(f" ({fmt(nm)} nm)" if nm is not None and unit!='nm' else '')
        pairs=[('Status',r['status']),('Accepted apparent Rg (nm)',fmt(r.get('rg_nm'))),
            ('Diagnostic candidate Rg; not validated if rejected',rgtext),('q fit interval',f"{fmt(c.get('qmin'))} to {fmt(c.get('qmax'))} {r['q_unit']}"),
            ('qmin Rg / qmax Rg',f"{fmt(c.get('qmin_rg'))} / {fmt(c.get('qmax_rg'))}"),('Fit points',fmt(c.get('n_points'))),
            ('R² (insufficient for validity)',fmt(c.get('r2'))),('Endpoint relative spread',fmt(c.get('endpoint_relative_spread'))),
            ('Reasons','; '.join(r.get('rejection_reasons',[])) or 'No automatic gate failed; physical assignment remains unconfirmed.'),
            ('Warnings','; '.join(r.get('warnings',[])))]
        table=''.join(f'<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>' for k,v in pairs)
        rep_sections.append(f'<details><summary>{esc(r["label"])} · {esc(r["status"])}</summary><table>{table}</table><img alt="Scattering and fit diagnostics" src="data:image/png;base64,{img}"></details>')
    src=''.join(f'<li><a href="{esc(u)}">{esc(t)}</a></li>' for t,u in SOURCES)
    narrative=esc((interpretation or {}).get('narrative','No external LLM was executed. Numerical results were calculated locally.'))
    content=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>PolymerScattLab · SAXS/WAXS report</title><style>
body{{font:16px/1.6 system-ui,sans-serif;background:#f2f5f8;color:#183043;margin:0}}main{{max-width:1100px;margin:auto;padding:40px 24px}}h1{{font-size:38px;line-height:1.2}}h2{{margin-top:34px}}.eyebrow{{color:#087e8b;font-weight:750;letter-spacing:.08em}}.cards{{display:flex;gap:14px;flex-wrap:wrap}}.card,section,details{{background:white;border:1px solid #dce4eb;border-radius:12px;padding:22px;margin:16px 0}}.card{{flex:1;min-width:150px}}.big{{font-size:32px;font-weight:750}}.notice{{border-left:5px solid #bd5a22;padding:16px 22px;background:#fff5e9}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;border-bottom:1px solid #e1e7ec;padding:10px;vertical-align:top}}th{{color:#486272}}img{{width:100%;height:auto}}a{{color:#087e8b}}summary{{font-weight:700;cursor:pointer}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font:14px/1.5 ui-monospace,monospace}}code{{background:#e9f1f4;padding:2px 5px}}.scroll{{overflow-x:auto}}@media print{{details{{display:block}}summary{{display:none}}main{{padding:0}}}}
</style><main><div class="eyebrow">POLYMERSCATTLAB · 0.1.0 · RESEARCH PROTOTYPE</div>
<h1>SAXS/WAXS analysis with auditable Rg screening</h1><p>Prepared for Bing Wu · {esc(manifest.get('created_utc',''))}</p>
<div class="cards"><div class="card"><div class="big">{counts['n_curves']:,}</div>profiles analyzed individually</div><div class="card"><div class="big">{counts['n_accepted']:,}</div>passing apparent Guinier screening</div><div class="card"><div class="big">{counts['n_rejected']+counts['n_invalid']:,}</div>rejected or unresolved</div></div>
<p class="notice"><strong>Molecular-chain Rg is not established.</strong> A numerical Guinier fit can describe a scattering feature, domain or aggregate. Polymer identity, physical state, concentration, contrast and background treatment are needed to assign it to individual chains. Rejected candidate values below are diagnostic slopes, not extracted polymer sizes.</p>
<section><h2>Sample results</h2><div class="scroll"><table><tr><th>Group</th><th>Curves</th><th>Pass</th><th>Reject</th><th>Invalid</th><th>q units</th></tr>{rows}</table></div><p>No frames were averaged. Counts indicate numerical screening outcomes, not distinct physical samples or independent replicates.</p>{notes}</section>
<section><h2>What the software does</h2><p>The numerical engine reads q, I and optional uncertainty, identifies no-data masks, screens a constrained set of low-q Guinier windows, calculates Rg and uncertainty, and checks range, statistical significance, residual curvature and endpoint stability. Every curve retains its status, reasons, source and diagnostic candidate. Files with unknown q units retain native lengths and cannot produce unconditional nanometre results.</p>
<p><code>ln I(q) = ln I(0) − q² Rg² / 3</code> and <code>Rg = √(−3 × slope)</code>. Supplied errors are used through σ(log I) ≈ σ(I)/I. Without errors, regression uncertainty measures scatter about the line only. It excludes calibration, subtraction, model selection, covariance between bins and physical-model uncertainty.</p>
<p>Defaults are conservative engineering checks: qmax·Rg ≤ 1, qmin·Rg ≤ 0.5, sufficient points and q-span, meaningful negative slope, no strong curvature and stable endpoints. These checks are heuristics, not a universal certification of Guinier validity. Only a small near-origin search is allowed to avoid finding accidental linear segments at high q. Candidate selection and thresholds are available in <code>core.py</code>.</p></section>
<section><h2>LLM interpretation layer</h2><pre>{narrative}</pre><p>The optional adapter can send compact fit summaries to a local Ollama model or an explicitly configured OpenAI-compatible endpoint. It explains results and requests missing metadata; it cannot modify the numerical result file. This is an inference integration, not an LLM trained on the uploaded experiments. Online provider execution is separate from numerical tests.</p></section>
<section><h2>WAXS and model extensions</h2><p>The package includes an exploratory peak/d-spacing function and an optional Debye ideal-coil fit. Peak spacing d = 2π/q is distinct from Rg. Peak detection requires a user-selected physical q interval and confirmation of units; it does not identify a crystalline phase, quantify crystallinity or estimate Scherrer sizes automatically. WAXS coverage must be checked from calibrated q limits and the material being studied; a filename alone does not establish it.</p>
<p>A Debye fit is conditional on ideal Gaussian-chain scattering and cannot rescue an unresolved Guinier result by itself. Constant background, form-factor amplitude and Rg can be highly correlated without a measured plateau. No forced Debye estimate is promoted to a molecular Rg in this report.</p></section>
<h2>Representative diagnostics</h2><p>First and last frame of each archive group, plus each workbook curve. All profile-level outputs are in the accompanying CSV/JSON.</p>{''.join(rep_sections)}
<section><h2>Information needed for physical interpretation</h2><ul><li>Polymer/sample mapping and whether each is a dilute solution, gel, film or bulk material.</li><li>Archive q units, instrument calibration, uncertainty-column definition and any desmearing.</li><li>Solvent, concentration, transmission normalization and sample/solvent/cell background protocol.</li><li>Meaning of frame indices and T0/T24/T48/T72; exposure duration and damage or processing history.</li><li>If possible, lower-q measurements and a dilution series to distinguish chains, interactions and aggregates.</li></ul></section>
<section><h2>Provenance and reproducibility</h2><pre>{esc(json.dumps(manifest,indent=2))}</pre><p>Input files are not modified. Reduced data are accepted as supplied; folder names suggest background subtraction but do not prove its quality.</p></section>
<section><h2>Primary and official method references</h2><ul>{src}</ul><p>References support the model definitions and interpretation. Acceptance thresholds are this prototype’s conservative design choices.</p></section></main></html>'''
    Path(path).write_text(content,encoding='utf-8')
