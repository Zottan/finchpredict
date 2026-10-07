# finchpredict
Genetics prediction for aviary finches

# AmanDNA | Díszmadár genetikai predikció – v 1.0

## Futtatás

Python 3.10+ szükséges, Python 3.12 ajánlott.

A mappában legyen együtt:

- `genetics_app.py`
- `mutations.json`

Futtatás:

```bash
python genetics_app.py
```

A program kizárólag a Python standard könyvtárát használja, ezért külön `pip install`
nem szükséges.

## Mutáció módosítása

A `mutations.json` szerkeszthető.

Példa:

```json
{
  "id": "opal",
  "name": "Opál",
  "symbol": "o",
  "inheritance": "autosomal_recessive",
  "lethal": false,
  "lethal_genotype": null
}
```

A `name` mezőben szabadon megváltoztatható a mutáció neve.

## Öröklésmenetek

- `autosomal_dominant`
- `autosomal_recessive`
- `z_linked_recessive`
- `z_linked_dominant`

Madár ivari rendszer:

- hím: ZZ
- tojó: ZW

## Következő fejlesztési irány

A későbbi ZebraCalc-szerű változatban érdemes különválasztani:

1. genetikai motor
2. fajadatbázis
3. mutációadatbázis
4. fenotípus-megjelenítés
5. madárrajz/renderelő rendszer
6. keresztezési fa
7. több gén egyidejű számítása
8. rekombináció/linkage
9. genetikai valószínűségek
10. menthető tenyésztési párok

A rajzokat nem feltétlenül kell egyenként megrajzolni. Célszerű lehet egy
réteges SVG-rendszer: test + fej + szárny + farok + szem + csőr + színréteg +
mintázat + mutációspecifikus réteg. Így egyetlen alapmadárból sok fenotípus
automatikusan előállítható.

