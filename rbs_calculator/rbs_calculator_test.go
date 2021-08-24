package rbs_calculator

var eColi16SrRNA string = Organism16SrRNAMap["Escherichia coli str. K-12 substr. DH10B"]

// `cd` into this directory and run with
// `go test -run ^ExampleRibosomeBindingSites | less -S`. Piping to
// `less -S` allows you to view the complete horizontal output (use the left
// and right arrow keys to scroll horizontally)
func ExampleRibosomeBindingSites() {
	startCodons := []StartCodon{AUG, GUG}
	mRNA := "UCUAGAGGCCGACGCAAGCCCAUAUCGGGGCUUCCGUCGGCCAUAAGGAGGUAAAAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAUGCGUUUCAAAGUUCGUAUGGAAGGUUCCGUUAACGGUCACGAGUUCGAAAUCGAAGGUGAAGGUGAA"
	ribosomeBindingSites := RibosomeBindingSites(eColi16SrRNA, mRNA, 37.0, startCodons)
	includeSequences, includeStructures := false, false
	PrintBindingSites(ribosomeBindingSites, includeSequences, includeStructures)

	// Output:
	// +----------------+-------------+---------------------+
	// | Start position | Start codon |         TIR         |
	// +----------------+-------------+---------------------+
	// |             57 | AUG         |   39800.46119720481 |
	// +----------------+-------------+---------------------+
	// |             90 | AUG         |   337.3785194542471 |
	// +----------------+-------------+---------------------+
	// |            108 | AUG         |   67.97945336767806 |
	// +----------------+-------------+---------------------+
	// |            148 | GUG         | 0.35638938610056775 |
	// +----------------+-------------+---------------------+
	// |            154 | GUG         |  10743.241317699827 |
	// +----------------+-------------+---------------------+
}

func ExampleRibosomeBindingSites_withSequences() {
	startCodons := []StartCodon{AUG, GUG}
	mRNA := "UCUAGAGGCCGACGCAAGCCCAUAUCGGGGCUUCCGUCGGCCAUAAGGAGGUAAAAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAUGCGUUUCAAAGUUCGUAUGGAAGGUUCCGUUAACGGUCACGAGUUCGAAAUCGAAGGUGAAGGUGAA"
	ribosomeBindingSites := RibosomeBindingSites(eColi16SrRNA, mRNA, 37.0, startCodons)
	includeSequences, includeStructures := true, false
	PrintBindingSites(ribosomeBindingSites, includeSequences, includeStructures)

	// Output:
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
	// | Start position | Start codon |         TIR         |                                        5' Untranslated Region                                        |       Protein Coding Sequence       |
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
	// |             57 | AUG         |   39800.46119720481 | UCUAGAGGCCGACGCAAGCCCAUAUCGGGGCUUCCGUCGGCCAUAAGGAGGUAAAAA                                            | AUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAU |
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
	// |             90 | AUG         |   337.3785194542471 | UCUAGAGGCCGACGCAAGCCCAUAUCGGGGCUUCCGUCGGCCAUAAGGAGGUAAAAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUC           | AUGCGUUUCAAAGUUCGUAUGGAAGGUUCCGUUAA |
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
	// |            108 | AUG         |   67.97945336767806 | CCGACGCAAGCCCAUAUCGGGGCUUCCGUCGGCCAUAAGGAGGUAAAAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAUGCGUUUCAAAGUUCGU | AUGGAAGGUUCCGUUAACGGUCACGAGUUCGAAAU |
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
	// |            148 | GUG         | 0.35638938610056775 | AGGUAAAAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAUGCGUUUCAAAGUUCGUAUGGAAGGUUCCGUUAACGGUCACGAGUUCGAAAUCGAAG | GUGAAGGUGAA                         |
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
	// |            154 | GUG         |  10743.241317699827 | AAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAUGCGUUUCAAAGUUCGUAUGGAAGGUUCCGUUAACGGUCACGAGUUCGAAAUCGAAGGUGAAG | GUGAA                               |
	// +----------------+-------------+---------------------+------------------------------------------------------------------------------------------------------+-------------------------------------+
}

func ExampleRibosomeBindingSites_withStructures() {
	startCodons := []StartCodon{AUG, GUG}
	mRNA := "UCUAGAGGCCGACGCAAGCCCAUAUCGGGGCUUCCGUCGGCCAUAAGGAGGUAAAAAAUGGCGAGCUCUGAAGACGUUAUCAAAGAGUUCAUGCGUUUCAAAGUUCGUAUGGAAGGUUCCGUUAACGGUCACGAGUUCGAAAUCGAAGGUGAAGGUGAA"
	ribosomeBindingSites := RibosomeBindingSites(eColi16SrRNA, mRNA, 37.0, startCodons)
	includeSequences, includeStructures := false, true
	PrintBindingSites(ribosomeBindingSites, includeSequences, includeStructures)

	// Output:
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
	// | Start position | Start codon |         TIR         |                                                              Initial state                                                              |                                      Final state (pre ribosome)                                      |    Final state (mRNA shine     | Final state (spacing) |     Final state (ribosome      | Final state (post ribosome) |  Final state (16S rRNA shine   |                                                             Full final state                                                             |
	// |                |             |                     |                                                                                                                                         |                                                                                                      |     dalgarno binding site)     |                       |           footprint)           |                             |     dalgarno binding site)     |                                                                                                                                          |
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
	// |             57 | AUG         |   39800.46119720481 | (((...((((((((.((((((......)))))).))))))))...)))..............(((((((...((.....))..)))))))..                                            | ......((((((((.((((((......)))))).)))))))).                                                          | .(((((((((                     | .....                 | .............                  | ......................      | &)))))))))                     | ......((((((((.((((((......)))))).))))))))..(((((((((........................................&)))))))))                                  |
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
	// |             90 | AUG         |   337.3785194542471 | (((...((((((((.((((((......)))))).))))))))...)))............(((((((.((((.((((...............)))))))).)))))))(((((....)))))...           | (((...((((((((.((((((......)))))).))))))))...)))..............(((((((...((.....))..)))))))           |                                |                       | .............                  | .....(((((....)))))...      | &                              | (((...((((((((.((((((......)))))).))))))))...)))..............(((((((...((.....))..)))))))..................(((((....)))))...&           |
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
	// |            108 | AUG         |   67.97945336767806 | ((((((.((((((......)))))).)))))).................((..((((((((((...(((((.(......((((((((..........)))))))).......).))))).))).)))))))..)) | ((((((.((((((......)))))).))))))....................(((((((.((((.((((...............)))))))).))))))) |                                |                       | .............                  | ......................      | &                              | ((((((.((((((......)))))).))))))....................(((((((.((((.((((...............)))))))).)))))))...................................& |
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
	// |            148 | GUG         | 0.35638938610056775 | ..........(..((((((((((...(((((.(......((((((((..........)))))))).......).))))).))).)))))))..)(((....))).......                         | ............(((((((.((((.((((...............)))))))).)))))))(((((....)))))....                       | .....(((                       | ...................   | .............                  |                             | &)))......                     | ............(((((((.((((.((((...............)))))))).)))))))(((((....))))).........(((................................&)))......         |
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
	// |            154 | GUG         |  10743.241317699827 | ....(..((((((((((...(((((.(......((((((((..........)))))))).......).))))).))).)))))))..)(((....))).......                               | ...((..((((((((((...(((((.(......((((((((..........)))))))).......).))))).))).)))))))..))...         | ......((((                     | ....                  | .............                  |                             | &)))).....                     | ...((..((((((((((...(((((.(......((((((((..........)))))))).......).))))).))).)))))))..)).........((((.................&)))).....        |
	// +----------------+-------------+---------------------+-----------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------+--------------------------------+-----------------------+--------------------------------+-----------------------------+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
}

// to run, cd into this directory and run with
// `go test -run ^ExampleRibosomeBindingSites_withPropertiesPrint | less -S`
func ExampleRibosomeBindingSites_withPropertiesPrint() {
	startCodons := []StartCodon{AUG, GUG, UUG}
	mRNA := "gctagcCACCGTCACACAGGAAAGtactagATGATTGAAAAAATTTGGAGCGGCGAAAGCCCGCATATGCGTAAAGGCGAAGAGCTGTTCACTGGTTTCGTCACTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCACATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCACATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATAA"
	ribosomeBindingSites := RibosomeBindingSites(eColi16SrRNA, mRNA, 37.0, startCodons)

	// make sure that each property included here is calculated for each RBS site
	// check the `PropertiesToCompute` variable (of your RBS calculator model
	// subpackage) to see the properties that can be printed
	propertiesToPrint := []RBSPropertyFunc{
		"dG_total", "dG_standby", "dG_mRNA_rRNA", "dG_spacing", "alignedSpacing", "dG_start", "dG_mRNA", "dG_stack", "spacingSequence"}
	includeSequences, includeStructures := false, false
	PrintBindingSites(ribosomeBindingSites, includeSequences, includeStructures, propertiesToPrint...)

	// Output:
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// | Start position | Start codon |         TIR          |      dG_total       |     dG_standby      |    dG_mRNA_rRNA     |      dG_spacing      | alignedSpacing | dG_start | dG_mRNA |       dG_stack        |   spacingSequence    |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |             30 | AUG         |   20183.660797413013 |   -6.21205151444689 |                   0 |  -10.99135151444689 |                    0 |              5 |    -2.76 |   -7.37 |                0.1693 | ACUAG                |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |             34 | UUG         |   1014.9135187674386 |  0.4325484855531103 |                   0 |  -10.99135151444689 |                1.728 |              9 |     1.81 |   -7.37 |                0.5159 | ACUAGAUGA            |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |             44 | UUG         |   1.0931721342382053 |   15.61804848555311 |                   0 | -7.8513515144468915 |               12.768 |             19 |     1.81 |   -7.46 |    1.4313999999999996 | ACUAGAUGAUUGAAAAAAU  |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |             66 | AUG         |    522.6194335885028 |  1.9074488098023785 |                   0 | -11.581351190197623 |                0.288 |              6 |    -2.76 |     -16 | -0.039199999999999985 | CGCAU                |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            111 | GUG         |  0.12295505895908687 |  20.473648595225654 |                   0 | -16.891351404774344 |    7.199999999999999 |             15 |    -0.42 |  -29.46 |    1.1249999999999998 | UUCGUCACUAUUCUG      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            121 | AUG         |    427.6622153416195 |  2.3530485237000818 |                   0 | -23.231351476299917 |                    0 |              5 |    -2.76 |  -28.13 |                0.2144 | AACUGG               |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            124 | GUG         |   18.215241114178067 |   9.366548523700082 |                   0 | -21.571351476299917 |   1.1520000000000001 |              8 |    -0.42 |  -29.69 |                0.5159 | AACUGGAUG            |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            127 | AUG         |    66.27141471043541 |   6.496548554694499 |                   0 |   -20.9913514453055 |                    0 |              5 |    -2.76 |  -29.69 |                0.5579 | AUGGUG               |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            150 | GUG         |    1.119305759255854 |  15.565548595225657 |                   0 |  -22.03135140477434 |                4.032 |             12 |    -0.42 |  -32.89 |                1.0949 | CAUAAGUUUUCC         |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            154 | GUG         |  0.07092364772816478 |  21.696348595225658 |                   0 | -17.771351404774343 |                8.448 |             16 |    -0.42 |   -30.2 |                1.2397 | CAUAAGUUUUCCGUGC     |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            163 | GUG         |   174.11317227505495 |   4.349999999999998 |                   0 |              -24.57 |                    0 |              0 |    -0.42 |  -29.34 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            169 | GUG         |   1025.1032829284668 | 0.41034856184705504 |                   0 | -28.901351438152943 |                    0 |              5 |    -0.42 |  -29.34 |   0.39170000000000005 | GUGAAG               |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            181 | AUG         |    602.9134373843755 |  1.5898485999940277 |                   0 |  -25.47135140000597 |                  2.4 |             10 |    -2.76 |  -27.43 | -0.008800000000000002 | GACGCAACUA           |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            233 | UUG         |   1.6740295685268582 |  14.671048595225656 |                   0 |  -19.53135140477434 |                1.525 |              3 |     1.81 |  -31.08 |               -0.2126 | ACC                  |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            262 | AUG         |   1.1024595466761558 |  15.599248595225657 |                   0 | -18.231351404774344 |                6.048 |             14 |    -2.76 |  -30.25 |   0.29259999999999997 | AACGACGCUGACUU       |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            265 | GUG         | 0.028048324969143876 |  23.757848595225656 |                   0 | -14.921351404774343 |    9.792000000000002 |             17 |    -0.42 |  -28.63 |    0.6771999999999999 | AACGACGCUGACUUAUG    |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            272 | GUG         |   3.9148954938590297 |  12.783148595225654 |                   0 | -15.721351404774344 |                    0 |              5 |    -0.42 |  -28.54 |   0.38449999999999995 | GUUCA                |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            277 | UUG         |   0.8797773941602054 |  16.100648595225657 |                   0 |  -17.58135140477434 |                  2.4 |             10 |     1.81 |  -28.79 |    0.6819999999999999 | GUUCAGUGCU           |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            297 | AUG         |    80.65699254430253 |   6.060000000000006 |                   0 | -22.549999999999997 |                    0 |              0 |    -2.76 |  -31.37 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            307 | AUG         |    9.462850397369388 |   10.82184856661543 |                   0 |  -23.59135143338457 |                3.168 |             11 |    -2.76 |  -34.06 |  -0.05480000000000003 | CCACAUGAAGCAGC       |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            327 | AUG         |    600.1655315960385 |  1.5999999999999979 |                   0 |              -27.48 |                    0 |              0 |    -2.76 |  -31.84 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            340 | AUG         |   19.131156993628405 |   9.257527202110584 |                   0 |  -29.30135119019762 |    9.628878392308204 |              2 |    -2.76 |  -31.69 |                     0 | U                    |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            342 | GUG         |    90.24579711678962 |   5.810374391943213 |                   0 | -30.461351190197625 | 0.005325582140839866 |              4 |    -0.42 |  -36.43 |                0.2564 | UAU                  |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            370 | AUG         |    0.536286710130751 |    17.2006485713838 |                   0 | -21.611351428616203 |    9.792000000000002 |             17 |    -2.76 |  -30.69 |    1.0899999999999996 | ACGCACGAUUUCCUUUAAGG |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            391 | GUG         |   1.3778422801410137 |   15.10374842833265 | 0.21680000000000277 | -23.021351571667353 |    9.792000000000002 |             17 |    -0.42 |  -28.54 | -0.003699999999999995 | UGACGGCACGUACAAAACGC |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            399 | GUG         |    696.2477106770248 |   1.269999999999996 |                   0 |              -26.01 |                    0 |              0 |    -0.42 |   -27.7 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            406 | UUG         |   359.19615771991147 |  2.7407485856889195 | 0.07860000000000511 | -26.941351414311086 |                    0 |              5 |     1.81 |  -27.53 |                0.2635 | GAAAU                |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            433 | UUG         |   25.478211484095986 |   8.620848595225656 |                   0 | -19.141351404774344 |   1.1520000000000001 |              8 |     1.81 |  -25.09 |               -0.2898 | AAACCGCA             |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            448 | UUG         |   1.9276822320577958 |  14.357526973228744 |                   0 |  -19.48135141907946 |    9.628878392308204 |              2 |     1.81 |   -22.4 |                     0 | A                    |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            511 | AUG         |   15.283955960721444 |   9.756448566615429 |                   0 | -17.101351433384572 |    9.792000000000002 |             17 |    -2.76 |     -19 |    0.8257999999999998 | AUACAAUUUUAACAGCCACA |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            541 | AUG         |   293.44948240033386 |  3.1899999999999995 |                   0 |               -9.72 |                    0 |              0 |    -2.76 |  -15.67 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            576 | GUG         |  0.12443567563986313 |   20.44704858062252 |                   0 | -12.981351419377482 |    9.792000000000002 |             17 |    -0.42 |  -22.89 |    1.1663999999999992 | AUUUUAAAAUUCGCCACAAC |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            583 | AUG         |    226.0385371638983 |  3.7699999999999996 |                   0 |              -17.98 |                    0 |              0 |    -2.76 |  -24.51 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            591 | GUG         |   223.02272588771177 |   3.799848504626592 |                   0 | -20.101351495373407 |   0.6719999999999999 |              7 |    -0.42 |  -23.47 |   0.17919999999999997 | AUGGCAGC             |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            631 | GUG         |   11.648827173716805 |               10.36 |                   0 |              -19.52 |                    0 |              0 |    -0.42 |   -30.3 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            634 | AUG         |      95.269098506299 |   5.690000000000001 |                   0 |              -21.85 |                    0 |              0 |    -2.76 |   -30.3 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            717 | AUG         |    398.4992049652895 |   2.510000000000005 |                   0 |              -23.65 |                    0 |              0 |    -2.76 |  -28.92 |                     0 |                      |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            757 | AUG         |    494.2115390759347 |   2.031648579728447 |                   0 | -19.691351420271552 |   0.6719999999999999 |              7 |    -2.76 |  -23.91 |                -0.099 | CAUCACGC             |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            762 | AUG         |    92.08468495331736 |   5.765548579728446 |                   0 |  -17.91135142027155 |                4.032 |             12 |    -2.76 |  -22.13 |   0.27490000000000003 | CAUCACGCAUGGU        |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
	// |            766 | AUG         |    2237.785687013095 | -1.3245258226335004 |                   0 | -20.971351404774342 | 0.005325582140839866 |              4 |    -2.76 |   -22.1 |                0.3015 | AUGG                 |
	// +----------------+-------------+----------------------+---------------------+---------------------+---------------------+----------------------+----------------+----------+---------+-----------------------+----------------------+
}
