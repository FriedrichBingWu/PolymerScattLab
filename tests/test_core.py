import numpy as np
import pytest
from polymer_scatt.core import analyze_curve,fit_debye,debye_form_factor,find_peaks


def test_guinier_recovers_size_and_native_units():
    q=np.linspace(.005,.14,100);I=30*np.exp(-(q*8)**2/3)
    r=analyze_curve(q,I,.005*I,q_unit="nm^-1")
    assert r["status"]=="accepted"
    assert r["rg_nm"]==pytest.approx(8,rel=1e-9)
    assert r["fit"]["qmax_rg"]<=1
    assert r["fit"]["qmin_rg"]<=.5
    assert analyze_curve(q,I)["rg_nm"] is None


def test_unit_scaling():
    q=np.linspace(.003,.12,90);I=40*np.exp(-(q*9)**2/3)
    a=analyze_curve(q,I,q_unit="nm^-1")
    b=analyze_curve(q/10,I,q_unit="A^-1")
    assert a["status"]==b["status"]=="accepted"
    assert a["rg_nm"]==pytest.approx(b["rg_nm"],rel=1e-9)
    assert b["rg_native"]==pytest.approx(a["rg_native"]*10,rel=1e-9)


def test_power_law_does_not_invent_plateau():
    q=np.geomspace(.02,2,200);I=q**-3
    r=analyze_curve(q,I,.02*I)
    assert r["status"]=="rejected"
    assert r["rg_native"] is None
    assert r["candidate"] is not None


def test_masks_and_invalid_errors():
    q=np.linspace(.001,.12,100);I=50*np.exp(-(q*7)**2/3);s=.01*I
    I[:8]=0;s[:8]=0
    r=analyze_curve(q,I,s)
    assert r["status"]=="accepted"
    assert r["diagnostics"]["zero_zero_masked_bins"]==8
    assert r["fit"]["start_index"]>=8
    s[25]=0
    assert analyze_curve(q,I,s)["rejection_reasons"]==["invalid_supplied_sigma"]


def test_no_positive_only_window_stitching():
    q=np.linspace(.005,.15,100);I=np.exp(-(q*7)**2/3)
    I[::4]=-1
    r=analyze_curve(q,I)
    assert r["status"]=="rejected"
    assert r["rg_native"] is None


def test_debye_recovery_and_stable_origin():
    assert debye_form_factor(np.array([0,1e-15]),8)==pytest.approx([1,1])
    q=np.geomspace(.004,1,160);I=23*debye_form_factor(q,8)
    r=fit_debye(q,I,.01*I,q_unit="nm^-1")
    assert r["status"]=="hypothesis_fit"
    assert r["rg_nm"]==pytest.approx(8,rel=1e-5)
    assert r["amplitude"]==pytest.approx(23,rel=1e-5)


def test_debye_background_requires_opt_in():
    q=np.geomspace(.004,1,160);I=23*debye_form_factor(q,8)+.4
    r=fit_debye(q,I,np.full(len(q),.01),q_unit="nm^-1",fit_background=True)
    assert r["rg_native"]==pytest.approx(8,rel=1e-4)
    assert r["background"]==pytest.approx(.4,rel=1e-4)


def test_debye_without_low_q_is_unresolved():
    q=np.geomspace(.5,3,100);I=23*debye_form_factor(q,8)
    r=fit_debye(q,I,.01*I)
    assert r["status"]=="unresolved"
    assert "low_q_plateau_not_resolved" in r["rejection_reasons"]


def test_internal_masks_do_not_join_disconnected_regions():
    q=np.linspace(.005,.15,100);I=np.exp(-(q*7)**2/3);s=.01*I
    I[8:11]=0;s[8:11]=0
    r=analyze_curve(q,I,s)
    assert r["status"]=="rejected"
    assert r["rg_native"] is None


def test_peaks_spacing_and_no_phase_claim():
    q=np.linspace(.1,3,300);I=4+20*np.exp(-((q-1.6)/.05)**2)
    r=find_peaks(q,I,q_unit="A^-1")
    assert len(r["peaks"])==1
    assert r["peaks"][0]["d_nm"]==pytest.approx(2*np.pi/1.6*.1,rel=.01)


def test_two_waxs_peaks_and_masked_detector_gap():
    q=np.linspace(.1,3,1201)
    I=3+20*np.exp(-((q-1.2)/.045)**2)+12*np.exp(-((q-2.1)/.06)**2)
    s=np.full(len(q),.1)
    gap=(q>1.55)&(q<1.65)
    I[gap]=0;s[gap]=0
    r=find_peaks(q,I,q_unit="A^-1",sigma=s)
    assert r["zero_zero_masked_bins"]==int(gap.sum())
    assert len(r["peaks"])==2
    for peak,position in zip(r["peaks"],[1.2,2.1]):
        assert peak["q_peak"]==pytest.approx(position,abs=.003)
        assert peak["d_native"]==pytest.approx(2*np.pi/position,rel=.003)
        assert peak["d_nm"]==pytest.approx(2*np.pi/position*.1,rel=.003)
        assert peak["half_prominence_width_q"]>0
        assert peak["index"] not in np.flatnonzero(gap)
    nm=find_peaks(q*10,I,q_unit="nm^-1",sigma=s)
    assert [p["d_nm"] for p in r["peaks"]]==pytest.approx([p["d_nm"] for p in nm["peaks"]])
    assert "crystallinity" not in r and "crystallite_size" not in r


def test_waxs_zero_masks_do_not_create_edge_peaks():
    q=np.linspace(.1,3,301);I=np.full(len(q),4.);s=np.ones(len(q))
    I[:25]=0;s[:25]=0;I[130:150]=0;s[130:150]=0
    r=find_peaks(q,I,q_unit="nm^-1",sigma=s)
    assert r["peaks"]==[]
    s[60]=0
    with pytest.raises(ValueError,match="Invalid supplied sigma"):
        find_peaks(q,I,q_unit="nm^-1",sigma=s)


def test_weak_flat_and_invalid_q():
    q=np.linspace(.01,1,100)
    assert analyze_curve(q,np.ones(100))["status"]=="rejected"
    assert analyze_curve(q[::-1],np.ones(100))["status"]=="invalid"
    with pytest.raises(ValueError):analyze_curve(q,np.ones(100),q_unit="guessed")
