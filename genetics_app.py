#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Díszmadár genetikai predikciós program
Python 3.10+ / 3.12
GUI: Tkinter
Adatbázis: mutations.json

Alap öröklésmenetek:
- autoszomális domináns
- autoszomális recesszív
- Z-kromoszómához kötött recesszív
- Z-kromoszómához kötött domináns
- letális faktorok

Madár ivari rendszer:
- hím: ZZ
- tojó: ZW

A program nem fajspecifikus. A mutations.json fájlban a mutációk
nevei és tulajdonságai szabadon módosíthatók.
"""

from __future__ import annotations

import json
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox


BASE_DIR = Path(__file__).resolve().parent
DB_FILE = BASE_DIR / "mutations.json"


@dataclass(frozen=True)
class Mutation:
    id: str
    name: str
    symbol: str
    inheritance: str
    lethal: bool = False
    lethal_genotype: str | None = None
    description: str = ""


def load_mutations() -> list[Mutation]:
    with DB_FILE.open("r", encoding="utf-8") as f:
        data = json.load(f)

    return [
        Mutation(
            id=m["id"],
            name=m["name"],
            symbol=m["symbol"],
            inheritance=m["inheritance"],
            lethal=m.get("lethal", False),
            lethal_genotype=m.get("lethal_genotype"),
            description=m.get("description", ""),
        )
        for m in data.get("mutations", [])
    ]


def weighted(items):
    """Összevonja az azonos eredményeket."""
    result = defaultdict(float)
    for genotype, probability in items:
        result[genotype] += probability
    return dict(result)


def gametes_autosomal(genotype):
    """
    Egy egyszerű autoszomális gén gamétái.
    Példák:
      AA -> A
      Aa -> A/a
      aa -> a
    """
    a, b = genotype
    if a == b:
        return {a: 1.0}
    return {a: 0.5, b: 0.5}


def sort_alleles(a, b):
    return "".join(sorted((a, b)))


def autosomal_cross(male_gt, female_gt):
    male = gametes_autosomal(male_gt)
    female = gametes_autosomal(female_gt)

    results = []
    for ma, mp in male.items():
        for fa, fp in female.items():
            results.append((sort_alleles(ma, fa), mp * fp))
    return weighted(results)


def z_recessive_cross(male_gt, female_gt):
    """
    Z-kapcsolt gén madaraknál.

    Hím: Z/Z
    Tojó: Z/W

    male_gt = tuple(Z allele, Z allele)
    female_gt = tuple(Z allele, W)

    Az utódokhoz explicit ivart is rendelünk.
    """
    male_gametes = gametes_autosomal(male_gt)
    female_gametes = gametes_autosomal((female_gt[0], "W"))

    results = []

    for male_z, mp in male_gametes.items():
        for female_chr, fp in female_gametes.items():
            p = mp * fp

            # Tojó: apa ad Z, anya ad W
            if female_chr == "W":
                gt = f"{male_z}/W"
                results.append(("FEMALE", gt, p))

            # Hím: anya ad Z, apa ad Z
            else:
                gt = sort_alleles(male_z, female_chr)
                results.append(("MALE", gt, p))

    result = defaultdict(float)
    for sex, gt, p in results:
        result[(sex, gt)] += p
    return dict(result)


def z_dominant_cross(male_gt, female_gt):
    return z_recessive_cross(male_gt, female_gt)


def is_lethal(mutation: Mutation, genotype: str, sex: str | None = None):
    if not mutation.lethal:
        return False

    if mutation.inheritance == "autosomal_dominant":
        if mutation.lethal_genotype == "homozygous":
            # csak akkor letális, ha mindkét allél mutáns
            return genotype.count(mutation.symbol) == 2

    if mutation.inheritance == "autosomal_recessive":
        if mutation.lethal_genotype in ("homozygous", "affected"):
            return genotype.count(mutation.symbol) == 2

    if mutation.inheritance.startswith("z_linked"):
        # Tojónál egyetlen mutáns Z is lehet érintett;
        # a konkrét letalitási modell később fajspecifikusan bővíthető.
        if sex == "FEMALE":
            return genotype.startswith(mutation.symbol + "/")
        return genotype.count(mutation.symbol) == 2

    return False


def phenotype(mutation: Mutation, genotype: str, sex: str | None = None):
    """
    Meghatározza a fenotípust az adott egyszerű modell szerint.
    """
    s = mutation.symbol

    if mutation.inheritance == "autosomal_recessive":
        return "mutáns" if genotype.count(s) == 2 else "normál / hordozó"

    if mutation.inheritance == "autosomal_dominant":
        return "mutáns" if s in genotype else "normál"

    if mutation.inheritance == "z_linked_recessive":
        if sex == "FEMALE":
            return "mutáns" if genotype.startswith(s + "/") else "normál"
        return "mutáns" if genotype.count(s) == 2 else "normál / hordozó"

    if mutation.inheritance == "z_linked_dominant":
        return "mutáns" if s in genotype else "normál"

    return "ismeretlen"


def make_genotype(symbol: str, state: str, sex: str, inheritance: str):
    """
    A GUI számára készít egyszerű genotípusokat.
    state:
      wild, carrier, affected, heterozygous
    """
    if inheritance.startswith("z_linked"):
        if sex == "MALE":
            if state == "wild":
                return ("+", "+")
            if state in ("carrier", "heterozygous"):
                return ("+", symbol)
            if state == "affected":
                return (symbol, symbol)
        else:
            if state == "wild":
                return ("+", "W")
            if state in ("carrier", "heterozygous"):
                # Z-kapcsolt recesszívnél ez már érintett állapot
                return (symbol, "W")
            if state == "affected":
                return (symbol, "W")

    else:
        if state == "wild":
            return ("+", "+")
        if state in ("carrier", "heterozygous"):
            return ("+", symbol)
        if state == "affected":
            return (symbol, symbol)

    raise ValueError("Érvénytelen genotípus-állapot.")


def calculate(mutation, male_state, female_state):
    inh = mutation.inheritance
    s = mutation.symbol

    if inh.startswith("z_linked"):
        male_gt = make_genotype(s, male_state, "MALE", inh)
        female_gt = make_genotype(s, female_state, "FEMALE", inh)
        raw = z_recessive_cross(male_gt, female_gt)

        results = []
        for (sex, gt), p in raw.items():
            lethal = is_lethal(mutation, gt, sex)
            results.append({
                "sex": "hím" if sex == "MALE" else "tojó",
                "genotype": gt,
                "phenotype": phenotype(mutation, gt, sex),
                "probability": p,
                "lethal": lethal,
            })
        return results

    male_gt = make_genotype(s, male_state, "MALE", inh)
    female_gt = make_genotype(s, female_state, "FEMALE", inh)
    raw = autosomal_cross(male_gt, female_gt)

    results = []
    for gt, p in raw.items():
        lethal = is_lethal(mutation, gt)
        results.append({
            "sex": "hím/tojó",
            "genotype": gt,
            "phenotype": phenotype(mutation, gt),
            "probability": p,
            "lethal": lethal,
        })

    # Autoszomális öröklésnél az ivar független: minden eredmény
    # fele hím, fele tojó. A GUI ezt külön mutatja.
    return results


class GeneticsApp(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Díszmadár genetikai predikció")
        self.geometry("950x650")
        self.minsize(850, 550)

        self.mutations = load_mutations()
        self.mutation_map = {m.name: m for m in self.mutations}

        self.create_widgets()
        self.update_mutation_info()

    def create_widgets(self):
        top = ttk.Frame(self, padding=12)
        top.pack(fill="x")

        ttk.Label(
            top,
            text="Díszmadár genetikai predikció",
            font=("TkDefaultFont", 18, "bold")
        ).pack(anchor="w")

        ttk.Label(
            top,
            text="Általános Mendel-modell – madaraknál ZZ hím / ZW tojó rendszerrel"
        ).pack(anchor="w", pady=(3, 12))

        form = ttk.LabelFrame(self, text="Keresztezés", padding=12)
        form.pack(fill="x", padx=12)

        ttk.Label(form, text="Mutáció:").grid(row=0, column=0, sticky="w")
        self.mutation_var = tk.StringVar()
        self.mutation_combo = ttk.Combobox(
            form,
            textvariable=self.mutation_var,
            values=list(self.mutation_map.keys()),
            state="readonly",
            width=35
        )
        self.mutation_combo.grid(row=0, column=1, padx=8, sticky="w")
        self.mutation_combo.bind("<<ComboboxSelected>>",
                                 lambda e: self.update_mutation_info())

        self.info_label = ttk.Label(form, text="", foreground="#444")
        self.info_label.grid(row=0, column=2, padx=15, sticky="w")

        ttk.Label(form, text="Hím:").grid(row=1, column=0, pady=(12, 0), sticky="w")
        self.male_var = tk.StringVar(value="wild")
        self.male_combo = ttk.Combobox(
            form,
            textvariable=self.male_var,
            values=["wild", "carrier", "heterozygous", "affected"],
            state="readonly",
            width=20
        )
        self.male_combo.grid(row=1, column=1, padx=8, pady=(12, 0), sticky="w")

        ttk.Label(form, text="Tojó:").grid(row=2, column=0, pady=8, sticky="w")
        self.female_var = tk.StringVar(value="wild")
        self.female_combo = ttk.Combobox(
            form,
            textvariable=self.female_var,
            values=["wild", "carrier", "heterozygous", "affected"],
            state="readonly",
            width=20
        )
        self.female_combo.grid(row=2, column=1, padx=8, pady=8, sticky="w")

        self.calculate_button = ttk.Button(
            form,
            text="KERESZTEZÉS KISZÁMÍTÁSA",
            command=self.run_calculation
        )
        self.calculate_button.grid(row=3, column=0, columnspan=2, pady=10, sticky="w")

        # Eredmény
        result_frame = ttk.LabelFrame(self, text="Várható utódok", padding=10)
        result_frame.pack(fill="both", expand=True, padx=12, pady=12)

        columns = ("sex", "genotype", "phenotype", "probability", "lethal")
        self.tree = ttk.Treeview(result_frame, columns=columns, show="headings")

        headings = {
            "sex": "Ivar",
            "genotype": "Genotípus",
            "phenotype": "Fenotípus",
            "probability": "Valószínűség",
            "lethal": "Letális"
        }

        widths = {
            "sex": 110,
            "genotype": 170,
            "phenotype": 220,
            "probability": 130,
            "lethal": 100
        }

        for col in columns:
            self.tree.heading(col, text=headings[col])
            self.tree.column(col, width=widths[col], anchor="center")

        scrollbar = ttk.Scrollbar(
            result_frame, orient="vertical", command=self.tree.yview
        )
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        self.summary = ttk.Label(
            self,
            text="",
            padding=(12, 0, 12, 12),
            justify="left"
        )
        self.summary.pack(fill="x")

    def update_mutation_info(self):
        mutation = self.mutation_map.get(self.mutation_var.get())
        if not mutation:
            return

        inheritance_names = {
            "autosomal_dominant": "autoszomális domináns",
            "autosomal_recessive": "autoszomális recesszív",
            "z_linked_recessive": "Z-kapcsolt recesszív",
            "z_linked_dominant": "Z-kapcsolt domináns",
        }

        text = (
            f"{inheritance_names.get(mutation.inheritance, mutation.inheritance)}"
            f" | jel: {mutation.symbol}"
        )

        if mutation.lethal:
            text += " | LETÁLIS FAKTOR"

        self.info_label.config(text=text)

    def run_calculation(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

        mutation = self.mutation_map.get(self.mutation_var.get())

        if not mutation:
            messagebox.showerror("Hiba", "Válassz mutációt.")
            return

        try:
            results = calculate(
                mutation,
                self.male_var.get(),
                self.female_var.get()
            )
        except Exception as e:
            messagebox.showerror("Genetikai számítási hiba", str(e))
            return

        total_live = 0
        total_lethal = 0

        for r in results:
            p = r["probability"] * 100

            if r["lethal"]:
                total_lethal += r["probability"]
                lethal_text = "IGEN"
            else:
                total_live += r["probability"]
                lethal_text = "nem"

            self.tree.insert(
                "",
                "end",
                values=(
                    r["sex"],
                    r["genotype"],
                    r["phenotype"],
                    f"{p:.2f} %",
                    lethal_text
                )
            )

        mutation_name = mutation.name

        if mutation.inheritance.startswith("z_linked"):
            note = (
                "Z-kapcsolt modell: hím = ZZ, tojó = ZW. "
                "A W kromoszóma nem hordozza a vizsgált Z-kapcsolt allélt."
            )
        else:
            note = (
                "Autoszomális modell: az ivar a gén öröklődésétől független, "
                "ezért az egyes genotípusokból várhatóan 50% hím és 50% tojó."
            )

        summary = (
            f"{mutation_name}: {mutation.inheritance}\n"
            f"Életképes utódok összesen: {total_live * 100:.2f} %\n"
            f"Letális kombinációk: {total_lethal * 100:.2f} %\n"
            f"{note}"
        )

        self.summary.config(text=summary)


if __name__ == "__main__":
    app = GeneticsApp()
    app.mainloop()
