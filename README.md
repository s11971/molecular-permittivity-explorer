# Molecular Permittivity Explorer

**Maxwell's Equations → Molecular Polarizability → Dielectric Constant**

무극성 분자의 polarizability와 분자 수밀도를 이용하여 Clausius–Mossotti 관계식으로 상대유전율 εᵣ을 계산하는 Streamlit 웹앱입니다.

## Local run

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Streamlit Community Cloud

GitHub 저장소를 선택하고 Main file path를 `app.py`로 지정하여 Deploy할 수 있습니다.

## 핵심 물리 모델

D = ε₀E + P

p = αE_loc

(εᵣ − 1)/(εᵣ + 2) = Nα/(3ε₀)

N = ρN_A/M

온도에 따른 밀도 ρ(T)를 이용하여 εᵣ(T)를 계산합니다.

## 주의

기체는 지정 압력에서 이상기체식으로 밀도를 계산하고, 액체는 기준 밀도와 체적 열팽창계수를 이용한 단순 1차 근사를 사용합니다. 실제 유전율은 주파수, 분자 간 상관, local-field approximation 등의 영향을 받으므로 이 앱은 실험값을 대체하기보다 분자 물성과 거시적 유전율의 이론적 연결을 탐구하는 모델입니다.
