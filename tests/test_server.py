"""In-process HTTP smoke tests; no model service or external internet required."""
import json
import re
import threading
from urllib import request,error
import pytest
from polymer_scatt.server import WorkbenchServer

@pytest.fixture
def server():
    srv=WorkbenchServer(('127.0.0.1',0))
    thread=threading.Thread(target=srv.serve_forever,daemon=True);thread.start()
    yield srv
    srv.shutdown();srv.server_close();thread.join()

def post(srv,path,body=None,headers=None):
    h={'Content-Type':'application/json','X-Workbench-Token':srv.token}
    if headers:h.update(headers)
    req=request.Request(f'http://127.0.0.1:{srv.server_port}{path}',
        data=json.dumps(body or {}).encode(),headers=h,method='POST')
    with request.urlopen(req) as response:return json.load(response)

def test_demo_fit_and_offline_explanation(server):
    loaded=post(server,'/api/demo')
    args={'session':loaded['session'],'curve':0,'q_unit':'nm^-1'}
    result=post(server,'/api/analyze',args)['result']
    assert result['status']=='accepted'
    assert abs(result['rg_nm']-5)<.05
    reply=post(server,'/api/explain',args)
    assert reply['executed'] is False
    assert 'no external or local language model was run' in reply['narrative']

def test_host_origin_token_and_file_isolation(server):
    with request.urlopen(f'http://127.0.0.1:{server.server_port}/') as response:
        text=response.read().decode()
    assert '__CSRF_TOKEN__' not in text and server.token in text
    for headers in [{'X-Workbench-Token':'bad'},{'Origin':'https://unrelated.example'},{'Host':'unrelated.example'}]:
        with pytest.raises(error.HTTPError) as exc:post(server,'/api/demo',headers=headers)
        assert exc.value.code==403
    with pytest.raises(error.HTTPError) as exc:
        request.urlopen(f'http://127.0.0.1:{server.server_port}/README.md')
    assert exc.value.code==404

def test_units_and_model_scope(server):
    loaded=post(server,'/api/demo');args={'session':loaded['session'],'curve':0}
    with pytest.raises(error.HTTPError) as exc:post(server,'/api/analyze',{**args,'q_unit':'A^-1'})
    assert exc.value.code==400
    post(server,'/api/analyze',args)
    with pytest.raises(error.HTTPError) as exc:
        post(server,'/api/explain',{**args,'provider':'ollama','endpoint':'http://external.example'})
    assert exc.value.code==400

def test_real_workbook_http_upload(server):
    from pathlib import Path
    path=Path(__file__).resolve().parents[1]/'examples'/'MAXS1.xlsx'
    req=request.Request(f'http://127.0.0.1:{server.server_port}/api/upload?name=MAXS1.xlsx&unit=unknown',
        data=path.read_bytes(),headers={'Content-Type':'application/octet-stream','X-Workbench-Token':server.token},method='POST')
    with request.urlopen(req) as response:loaded=json.load(response)
    assert [x['label'] for x in loaded['curves']]==['T0','T24','T48','T72']
    for idx in range(4):
        r=post(server,'/api/analyze',{'session':loaded['session'],'curve':idx})['result']
        assert r['status']=='rejected' and r['rg_nm'] is None
