# Díszmadár genetikai kalkulátor

Kísérleti, fajfüggetlen genetikai motor és fenotípus-renderelő.

## Részek
- `genetics_app.py` – Python/Tkinter prototípus
- `mutations_schema_extended.json` – bővített mutációs adatmodell
- `web/` – GitHub Pages-kompatibilis böngészős UAT
- `web/amandina_uat.svg` – rétegezett ékfarkú amandina tesztmodell

## Böngészős UAT
A `web/index.html` statikus JavaScript alkalmazás. A test, has, fej, szárny, szárnyminta, farok, csőr és láb színe élőben módosítható, a módosított SVG exportálható.

## Tervezett architektúra
genotípus → genetikai motor → fenotípus → vizuális paraméterek → SVG

A következő nagy lépés a több-lokuszos genetikai motor és a webes UAT összekötése.
