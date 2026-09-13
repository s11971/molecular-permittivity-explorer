import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Nonpolar Molecule -> Relative Permittivity Calculator
# Model:
#   Clausius-Mossotti:
#   (eps_r - 1)/(eps_r + 2) = N * alpha_SI / (3 eps0)
#
# alpha is stored as molecular polarizability volume in A^3.
# Conversion:
#   alpha_SI = 4*pi*eps0*alpha_A3*1e-30
#
# For gases:
#   rho(T,P) = P*M/(R*T)          [ideal gas]
#
# For liquids in this MVP:
#   rho(T) = rho_ref * [1 - beta*(T-T_ref)]
#   This is intentionally a simple approximation, not a full
#   equation-of-state or tabulated NIST density model.
# ------------------------------------------------------------

st.set_page_config(
    page_title="Nonpolar Molecule Permittivity",
    page_icon="⚡",
    layout="wide"
)

EPS0 = 8.8541878128e-12
NA = 6.02214076e23
R = 8.31446261815324

DATA_FILE = Path(__file__).with_name("nonpolar_molecules.json")
data = json.loads(DATA_FILE.read_text(encoding="utf-8"))
df = pd.DataFrame(data)

st.title("⚡ 무극성 분자의 유전율 계산기")
st.caption("Maxwell 물질방정식 → 분극 → Clausius–Mossotti → 분자 물성")

with st.expander("이 앱이 실제로 계산하는 물리 모델", expanded=False):
    st.latex(r"\mathbf D=\epsilon_0\mathbf E+\mathbf P")
    st.latex(r"\mathbf p=\alpha\mathbf E_{\rm loc}")
    st.latex(r"\frac{\epsilon_r-1}{\epsilon_r+2}=\frac{N\alpha}{3\epsilon_0}")
    st.markdown("""
**핵심:** 무극성 분자는 영구 쌍극자 모멘트가 0이지만 외부 전기장에 의해
유도 쌍극자를 만들 수 있습니다. 이 앱에서는 그 반응을 분자
polarizability \\(\\alpha\\)로 나타내고, 분자 수밀도 \\(N\\)을 곱하여
거시적 분극과 상대유전율을 연결합니다.

온도는 무극성 분자에 대한 이 단순 모델에서 주로 **밀도 변화**를 통해
들어갑니다.
""")

st.sidebar.header("계산 조건")
query = st.sidebar.text_input("분자 검색", placeholder="예: benzene, CH4, 71-43-2")
filtered = df.copy()
if query.strip():
    q = query.strip().lower()
    mask = (
        filtered["name"].str.lower().str.contains(q, regex=False)
        | filtered["formula"].str.lower().str.contains(q, regex=False)
        | filtered["cas"].str.lower().str.contains(q, regex=False)
    )
    filtered = filtered[mask]

if filtered.empty:
    st.error("검색 결과가 없습니다. 이름, 분자식, CAS RN 중 하나로 검색해 보세요.")
    st.stop()

labels = [
    f"{row['name']} ({row['formula']}) — CAS {row['cas']}"
    for _, row in filtered.iterrows()
]
selected_label = st.sidebar.selectbox("분자 선택", labels)
selected = filtered.iloc[labels.index(selected_label)]

st.sidebar.divider()
pressure_bar = st.sidebar.number_input(
    "압력 (bar) — 기체 계산용", min_value=0.01, max_value=100.0,
    value=1.01325, step=0.1
)

tmin, tmax = st.sidebar.slider(
    "온도 범위 (°C)", min_value=-100, max_value=300,
    value=(-20, 100), step=5
)
npoints = st.sidebar.slider("그래프 계산점", 20, 300, 120)

def density_g_cm3(row, temp_C, pressure_bar):
    T = np.asarray(temp_C) + 273.15
    if row["phase"] == "gas":
        P = pressure_bar * 1e5
        M_kg_mol = row["molar_mass_g_mol"] / 1000.0
        rho_kg_m3 = P * M_kg_mol / (R * T)
        return rho_kg_m3 / 1000.0
    else:
        rho0 = row["density_ref_g_cm3"]
        T0 = row["density_ref_C"]
        beta = row["thermal_expansion_per_K"]
        rho = rho0 * (1.0 - beta * (np.asarray(temp_C) - T0))
        return np.maximum(rho, 1e-6)

def permittivity_from_density(row, rho_g_cm3):
    rho_kg_m3 = np.asarray(rho_g_cm3) * 1000.0
    M_kg_mol = row["molar_mass_g_mol"] / 1000.0
    N = rho_kg_m3 / M_kg_mol * NA

    alpha_A3 = row["polarizability_A3"]
    alpha_SI = 4 * math.pi * EPS0 * alpha_A3 * 1e-30

    A = N * alpha_SI / (3 * EPS0)
    if np.any(A >= 1):
        return np.full_like(np.asarray(A, dtype=float), np.nan)

    return (1 + 2*A) / (1 - A)

temps = np.linspace(tmin, tmax, npoints)
rho = density_g_cm3(selected, temps, pressure_bar)
eps = permittivity_from_density(selected, rho)

T_display = 25.0 if tmin <= 25 <= tmax else (tmin + tmax) / 2
rho_display = float(density_g_cm3(selected, T_display, pressure_bar))
eps_display = float(permittivity_from_density(selected, rho_display))

c1, c2, c3, c4 = st.columns(4)
c1.metric("분자", selected["formula"])
c2.metric("Polarizability", f"{selected['polarizability_A3']:.3g} Å³")
c3.metric(f"밀도 @ {T_display:.1f} °C", f"{rho_display:.5g} g/cm³")
c4.metric(f"εᵣ @ {T_display:.1f} °C", f"{eps_display:.9f}")

st.write(
    f"**{selected['name']}** · {selected['formula']} · CAS RN {selected['cas']} · "
    f"상(25 °C 기준 데이터셋): **{selected['phase']}**"
)

left, right = st.columns(2)
with left:
    fig1, ax1 = plt.subplots(figsize=(7, 4))
    ax1.plot(temps, rho, linewidth=2)
    ax1.set_xlabel("Temperature (°C)")
    ax1.set_ylabel("Density (g/cm³)")
    ax1.set_title("Temperature–Density")
    ax1.grid(True, alpha=0.25)
    st.pyplot(fig1, clear_figure=True)

with right:
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    ax2.plot(temps, eps, linewidth=2)
    ax2.set_xlabel("Temperature (°C)")
    ax2.set_ylabel("Relative permittivity εᵣ")
    ax2.set_title("Temperature–Permittivity (Clausius–Mossotti)")
    ax2.grid(True, alpha=0.25)
    st.pyplot(fig2, clear_figure=True)

with st.expander("계산 과정 보기", expanded=True):
    rho25 = float(density_g_cm3(selected, 25.0, pressure_bar))
    M = selected["molar_mass_g_mol"]
    alpha_A3 = selected["polarizability_A3"]
    alpha_SI = 4 * math.pi * EPS0 * alpha_A3 * 1e-30
    N25 = rho25 * 1000 / (M / 1000) * NA
    A25 = N25 * alpha_SI / (3 * EPS0)
    eps25 = (1 + 2*A25) / (1 - A25)

    st.latex(r"\alpha_{\rm SI}=4\pi\epsilon_0\alpha_{\AA^3}\times10^{-30}")
    st.latex(r"N=\frac{\rho N_A}{M}")
    st.latex(r"\frac{\epsilon_r-1}{\epsilon_r+2}=\frac{N\alpha}{3\epsilon_0}")
    st.write(f"- 25 °C 밀도: `{rho25:.8g} g/cm³`")
    st.write(f"- 분자 수밀도 N: `{N25:.6e} m⁻³`")
    st.write(f"- α (SI): `{alpha_SI:.6e} F·m²`")
    st.write(f"- A = Nα/(3ε₀): `{A25:.6e}`")
    st.write(f"- 계산된 εᵣ: **{eps25:.9f}**")

st.divider()
st.subheader("여러 무극성 분자 비교")
compare_labels = [f"{row['name']} ({row['formula']})" for _, row in df.iterrows()]
compare = st.multiselect(
    "비교할 분자 선택 (최대 6개)",
    compare_labels,
    default=[f"{selected['name']} ({selected['formula']})"]
)
compare = compare[:6]

if compare:
    fig3, ax3 = plt.subplots(figsize=(10, 5))
    for label in compare:
        row = df[(df["name"] + " (" + df["formula"] + ")") == label].iloc[0]
        rho_i = density_g_cm3(row, temps, pressure_bar)
        eps_i = permittivity_from_density(row, rho_i)
        ax3.plot(temps, eps_i, linewidth=2, label=row["formula"])
    ax3.set_xlabel("Temperature (°C)")
    ax3.set_ylabel("Relative permittivity εᵣ")
    ax3.set_title("Nonpolar Molecule Comparison")
    ax3.grid(True, alpha=0.25)
    ax3.legend()
    st.pyplot(fig3, clear_figure=True)

with st.expander("내장 데이터셋 30종"):
    display_df = df[[
        "name", "formula", "cas", "phase",
        "molar_mass_g_mol", "polarizability_A3",
        "density_ref_g_cm3", "density_ref_C", "thermal_expansion_per_K"
    ]].copy()
    display_df.columns = [
        "Name", "Formula", "CAS RN", "Phase", "M (g/mol)",
        "α (Å³)", "ρ_ref (g/cm³)", "T_ref (°C)", "β (1/K)"
    ]
    st.dataframe(display_df, width='stretch', hide_index=True)

st.info(
    "주의: 이 버전은 연구용 MVP입니다. 기체는 이상기체식으로 밀도를 계산하고, "
    "액체는 기준 밀도와 체적 열팽창계수를 이용한 1차 근사식을 사용합니다. "
    "따라서 실제 측정 유전율과의 비교에서는 분자 간 상관, 분산력, 주파수 의존성, "
    "Clausius–Mossotti 근사의 한계를 반드시 논의해야 합니다."
)

st.caption(
    "내장 polarizability 값은 대표적인 문헌/데이터베이스 값을 사용한 계산용 데이터셋이며, "
    "정밀한 연구에서는 각 물질·온도·상태에 맞는 원자료를 확인해 교체하는 것을 권장합니다."
)
