"""Collection of experimentally derived rates for n-alkyl radicals"""

import numpy as np
import sympy

h_planck = 4.135667696e-15  # eV/s
R_kcal = 1.987e-3  # Gas constant in kcal/(mol*K)
R_cal = 1.987
exp = np.exp
ln = np.log


class Paper:
    def __init__(self,
                 authors="",
                 year="",
                 title="",
                 type="",
                 link="",
                 id=-1,
                 ):
        self.authors = authors
        self.link = link
        self.year = year
        self.title = title
        self.type = type
        self.id = id


class RateResult:
    def __init__(self,
                 rate_expression,
                 paper_id,
                 temp_range=[np.nan, np.nan],
                 reaction_type="1-x",
                 base_molecule="alkyl",
                 E_unit="n/a",
                 ):
        assert "=" in rate_expression
        k_part, rate_equation = rate_expression.split("=")
        self.rate_expression = rate_equation
        self.paper_id = paper_id
        self.reaction_type = reaction_type
        self.base_molecule = base_molecule
        self.E_unit = E_unit
        self.temp_range = temp_range

        assert self.reaction_type != "1-x"

        self.k_part = k_part

        if k_part == "k":
            self.basis = lambda x: x
        elif k_part == "log(k)":
            self.basis = lambda x: 10**x
        elif k_part == "ln(k)":
            self.basis = lambda x: np.exp(x)
        else:
            raise ValueError(f"Unrecognized k_part {k_part}")

    def __str__(self):
        return f"ID:{self.paper_id} ({self.base_molecule}: {self.reaction_type}, {self.temp_range[0]}-{self.temp_range[1]}K)"

    def latex(self):
        T = sympy.Symbol('T')
        rate_at_500K = self.get_rate(500)
        return rf"{self.reaction_type} & {self.base_molecule} & ${self.k_part}={sympy.latex(sympy.sympify(self.rate_expression))}$ & \num{{{rate_at_500K:.1e}}} & {self.temp_range[0]}-{self.temp_range[1]} & {{\tiny ID:{self.paper_id}}} \\"

    def get_rate(self, T):
        if self.E_unit == "cal":
            R = R_cal
        elif self.E_unit == "kcal":
            R = R_kcal
        elif self.E_unit == "n/a":
            R = None
        else:
            raise ValueError(f"Unrecognized E_unit {self.E_unit}")
        res = eval(self.rate_expression)
        res = self.basis(res)
        return res


rate_results = [
    RateResult(
        paper_id=5,
        base_molecule="hexyl",
        reaction_type="1-5",
        rate_expression="log(k)=9.41-(11.2e3/(2.3*R*T))",
        temp_range=[25+273.15, 105+273.15],
        E_unit="cal"
    ),
    RateResult(
        paper_id=6,
        base_molecule="hexyl",
        reaction_type="1-5",
        rate_expression="ln(k)=9.5*ln(10)-11.6/(R*T)",
        temp_range=[300,453],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=7,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="k=4.88e8*(T**0.846)*exp(-19.53/(R*T))",
        temp_range=[350,1300],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=7,
        base_molecule="hexyl",
        reaction_type="1-5",
        rate_expression="k=6.65e7*(T**0.823)*exp(-12.45/(R*T))",
        temp_range=[350,1300],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=8,
        base_molecule="hexyl",
        reaction_type="1-5",
        rate_expression="k=1.83e2*(T**2.55)*exp(-5516/T)",
        temp_range=[500, 1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=8,
        base_molecule="hexyl",
        reaction_type="1-4",
        rate_expression="k=6.98*(T**3.2)*exp(-8333/T)",
        temp_range=[500, 1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=9,
        base_molecule="heptyl",
        reaction_type="1-5",
        rate_expression="k=10**(2.83)*(T**2.39)*exp(-5237/T)",
        temp_range=[860,1050],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=9,
        base_molecule="heptyl",
        reaction_type="1-4",
        rate_expression="k=10**(2.07)*(T**2.85)*exp(-8680/T)",
        temp_range=[860,1050],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=9,
        base_molecule="heptyl",
        reaction_type="1-6",
        rate_expression="k=10**(2.39)*(T**2.51)*exp(-6292/T)",
        temp_range=[860,1050],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=13,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="k=3.3e8*exp(-15.1/(R*T))",
        temp_range=[24+273.15, 162+273.15],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=14,
        base_molecule="hexyl",
        reaction_type="1-4",
        rate_expression="log(k)=11.0-(24/(2.3*R*T))",
        temp_range=[723,823],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=14,
        base_molecule="hexyl",
        reaction_type="1-5",
        rate_expression="log(k)=10.5-(17/(2.3*R*T))",
        temp_range=[723,823],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=17,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="k=1.4e7*exp(-10.8e3/(R*T))",
        temp_range=[438.5,502.5],
        E_unit="cal"
    ),
    RateResult(
        paper_id=18,
        base_molecule="hexyl",
        reaction_type="1-5",
        rate_expression="k=2.0e7*exp(-8.3e3/(R*T))",
        temp_range=[350,410],
        E_unit="cal"
    ),
    RateResult(
        paper_id=19,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="log(k)=11.08-20.04/(2.303*R*T)",
        temp_range=[438,923],
        E_unit="kcal"
    ),
    RateResult(
        paper_id=20,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="k=2.432e3*(T**2.324)*exp(-8183/T)+9.111e5*exp(-5303/T)",
        temp_range=[300,1300],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=22,
        base_molecule="octyl",
        reaction_type="1-4",
        rate_expression="k=10**(0.71)*(T**3.23)*exp(-8479/T)",
        temp_range=[700,1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=22,
        base_molecule="octyl",
        reaction_type="1-5",
        rate_expression="k=10**(1.36)*(T**2.82)*exp(-5413/T)",
        temp_range=[700,1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=22,
        base_molecule="octyl",
        reaction_type="1-6",
        rate_expression="k=10**(0.47)*(T**3.08)*exp(-5544/T)",
        temp_range=[700,1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=24,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="k=10**(8.57)*((T/298)**3.030)*exp(-7696/T)",
        temp_range=[700,1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=24,
        base_molecule="pentyl",
        reaction_type="1-3",
        rate_expression="k=10**(5.53)*((T/298)**6.837)*exp(-9444/T)",
        temp_range=[700,1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=25,
        base_molecule="pentyl",
        reaction_type="1-4",
        rate_expression="k=10**(1.06)*(T**3.033)*exp(-7706/T)",
        temp_range=[400,1900],
        E_unit="n/a"
    ),
    RateResult(
        paper_id=25,
        base_molecule="pentyl",
        reaction_type="1-3",
        rate_expression="k=10**(-11.41)*(T**6.843)*exp(-9451/T)",
        temp_range=[400,1900],
        E_unit="n/a"
    ),
]
for r in rate_results:
    print(f"{r.get_rate(500):.2e} {str(r)}")


results_by_shift = lambda shft, rs: [r for r in rs if r.reaction_type==shft]


rate_results_all = {
    "1-3": None,
    "1-4": None,
    "1-5": None,
    "1-6": None,
}

rate_results_all = {shft: results_by_shift(shft, rate_results) for shft in rate_results_all.keys()}

