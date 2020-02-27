root -b -q smallntuple_chain.C\(\"DY\"\)
root -b -q smallntuple_chain.C\(\"TT\"\)
root -b -q smallntuple_chain.C\(\"WW\"\)
root -b -q smallntuple_chain.C\(\"ST\"\)
root -b -q smallntuple_chain.C\(\"ZP2000\"\)
root -b -q smallntuple_chain.C\(\"ZP2500\"\)
root -b -q smallntuple_chain.C\(\"ZP3000\"\)
hadd zprimetraining.root small*.root
rm small*.root

