package salis_lab_v2_1

import (
	"math"
	"sort"

	"github.com/TimothyStiles/poly/energy_params"
	"github.com/TimothyStiles/poly/linearfold"
	"github.com/TimothyStiles/poly/mfe"
	. "github.com/TimothyStiles/poly/rbs_calculator/model"
	hm "github.com/TimothyStiles/poly/rbs_calculator/model/salis_lab_v2_1/hairpin_modules"
	"github.com/TimothyStiles/poly/rbs_calculator/model/salis_lab_v2_1/shine_dalgarno_binding_site"
	. "github.com/TimothyStiles/poly/secondary_structure"
)

/******************************************************************************

The Salis Lab v2.1 RBS Calculator Model

This model uses a free energy model which quantifies the thermodynamic
interactions in terms of Gibbs free energy changes. For more details of the
model, please have a look at:
* https://docs.denovodna.com/docs/rbs-calculator/rbs-calculator-free-energy-model
* https://github.com/barricklab/ostir/wiki/Background (This explains the Salis
	Lab v1 RBS Calculator Model. Note that the standby free energy terms is
	computed differently in v2.1 compared to v1, and the constraints are
	different. For a detailed explanation of the differences between the v1 and
	v2.1 models, please have a look at the `Supplementary Tables, Figures, and
	Discussion` pdf at https://pubs.acs.org/doi/10.1021/acssynbio.0c00394
	(Supplementary Discussion 1 on page 12).

******************************************************************************/

/***********************
// Steps:
// 1. Set mRNA to 5p UTR + first 35 nts of CDS
// -> Fold with RNAfold (no dangling and 2007 params) to get dG_mRNA
// 2. Cofold 5p UTR with 16S rRNA
// -> Note down most 5p and 3p of binding site
// 3. Get dG_pre_ribosome by folding from start of 5pUTR to 5p annd then using
// RNAeval to calculate dG.
// 4. Calculate dG_SD_aSD by passing only basepairs folded to aSD to RNAeval. Add
// _calc_dG_init_adjustment to the value.
// 5. Fold the post ribosome sequences (nts 14 to 35 of CDS) to get dG_post_ribosome
// 6. Add dG_start based on start codon lookup value.
// 7. Add dG_spacing based on aligned spacing between SD and CDS.
// 8. Pass pre_ribosome sequence to get all hairpin modules.
// 9. For each module, successively unfold pairs from the base of the hairpins.
// Keep track of the energy require to break that base pair as dG_unfolding and
// keep track of the dG_distortion (energy due to amount of available surface area).
// Choose dG_distortion and dG_unfolding that minimises the sum of the two.
// 10. For each module, get the height of all the downstream hairpins and add a
// penalising dG_sliding term to the module.
// 11. Choose the module that minimises the dG_standby term.
// 12. Compute dG_total and translation initaltion rate


Functions to compute properties of mRNAs
***********************/

// PropertiesToCompute list the properties that need to be computed for a
// `model.RibosomeBindingSite` struct to calculate its translation initiation
// rate.
var PropertiesToCompute []RBSPropertyFunc = []RBSPropertyFunc{
	ensureSaneInput,

	// Step 1: Initialization of `rbs` struct
	StartPosition,
	FivePrimeUTRCutoff,
	updateFivePrimeUTR,
	CodingSequenceCutoff,
	updateCodingSequence,

	// Step 2: Calculate dG_start
	CDSStartCodonFreeEnergy,

	// Step 3: Find the Shine-Dalgarno Binding Site
	shineDalgarnoBindingSite,
	HasShineDalgarnoBindingSite,
	SDBindingSiteMRNAStructure,
	SDBindingSiteRRNAStructure,
	BindingSiteMRNAStructure,
	LenSDSequence,
	BindingSiteFivePrimeIdx,
	BindingSiteThreePrimeIdx,

	// Step 4: Calculate dG_spacing
	SpacingRegionFreeEnergy,
	AlignedSpacing,

	// Step 5: Calculate the free energy of the mRNA sequence (dG_mRNA)
	UsedMRNASequence,
	UsedMRNADotBracketStructure,
	MRNAFreeEnergy,

	// Step 6: Calculate dG_mRNA_rRNA (dG_pre_ribosome + dG_SD_aSD + dG_post_ribosome)
	// Calculate dG_pre_ribosome
	MRNAPreRibosomeSequence,
	MRNAPreRibosomeDotBracketStructure,
	PreRibosomeFreeEnergy,

	// Calculate dG_SD_aSD
	ShineDalgarnoHybridizationFreeEnergy,

	// Calculate dG_post_ribosome
	MRNAPostRibosomeSequence,
	MRNAPostRibosomeDotBracketStructure,
	PostRibosomeFreeEnergy,

	MRNARRNAHybridizationFreeEnergy,

	// Step 7: Calculate dG_standby
	MRNAPreRibosomeSecondaryStructure,
	HairpinModule,
	HasHairpinModule,
	StandbyModuleTotalFreeEnergy,
	// The rest is only for output
	StandbyModuleUnfoldingFreeEnergy,
	StandbyModuleSlidingFreeEnergy,
	StandbyModuleDistortionFreeEnergy,
	HairpinModuleDistalLength,
	HairpinModuleProximalLength,
	HairpinModuleHairpinHeight,
	HairpinModuleSurfaceArea,
	HairpinModuleFivePrimeIdx,
	HairpinModuleThreePrimeIdx,

	// Step 8: Calculate dG_stack
	SpacingSequence,
	LenSpacingSequence,
	SpacingSequenceStackFreeEnergy,

	TotalFreeEnergy,
}

// ShineDalgarnoBindingSite is a type alias for
// `shine_dalgarno_binding_site.ShineDalgarnoBindingSite`
type ShineDalgarnoBindingSite = shine_dalgarno_binding_site.ShineDalgarnoBindingSite

// read the docs of `shine_dalgarno_binding_site.LookupTable()` for more information
var sdBindingSiteLookupTable map[string]map[string]map[float64]ShineDalgarnoBindingSite

// init is run automagically whenever this package is called.
// https://www.digitalocean.com/community/tutorials/understanding-init-in-go
func init() {
	sdBindingSiteLookupTable = shine_dalgarno_binding_site.LookupTable()
}

/*******************************************************************************
Start of properties
*******************************************************************************/

/*****************************************************
Step 1: Initialization of `rbs` struct

This section sets the cutoff values for the
`FivePrimeUTR` and `ProteinCodingSequence` fields of
the `rbs` struct and updates the fields based on their
cutoff values
*****************************************************/

// ensureSaneInput ensures sane input
func ensureSaneInput(rbs *RibosomeBindingSite) interface{} {
	if len(rbs.ProteinCodingSequence) == 0 || len(rbs.RibosomalRNA) == 0 {
		panic("the `RibosomeBindingSite` has either no `ProteinCodingSequence` or `RibosomalRNA` set which is required to calculate translation initiation rate")
	}
	return nil
}

// StartPosition makes note of the position at which the protein coding sequence
// starts in the mRNA sequence. This property is not used in the free energy
// model, but is added as we need to output the position at which the mRNA
// sequence is delimited into its five prime UTR and coding sequence.
func StartPosition(rbs *RibosomeBindingSite) interface{} {
	return len(rbs.FivePrimeUTR)
}

// FivePrimeUTRCutoff sets the cutoff of the five prime UTR (from the 3' to 5'
// direction) to 100 (as specified by The Salis Lab RBS Calc v2.1 model) or to
// the length of the five prime UTR if its length is less than 100.
func FivePrimeUTRCutoff(rbs *RibosomeBindingSite) interface{} {
	fivePrimeUTRCutoff := 100

	lenFivePrimeUTR := len(rbs.FivePrimeUTR)
	if fivePrimeUTRCutoff > lenFivePrimeUTR {
		fivePrimeUTRCutoff = lenFivePrimeUTR
	}
	return fivePrimeUTRCutoff
}

// updateFivePrimeUTR cuts off the five prime UTR from the 3' to 5'
// direction to the value of the `fivePrimeUTRCutoff` property value, and
// updates the `FivePrimeUTR` field of `rbs`. NOTE: as principle, properties
// shouldn't update the fields of the `rbs` struct, but this is an
// initialization step so is an exception.
func updateFivePrimeUTR(rbs *RibosomeBindingSite) interface{} {
	fivePrimeUTRCutoff := rbs.PropertyValue(FivePrimeUTRCutoff).(int)
	rbs.FivePrimeUTR = rbs.FivePrimeUTR[len(rbs.FivePrimeUTR)-fivePrimeUTRCutoff:]
	return nil
}

// CodingSequenceCutoff sets the cutoff of the coding sequence (from the 5' to
// 3' direction) to 35 (as specified by The Salis Lab RBS Calc v2.1 model) or to
// the length of the coding sequence if its length is less than 35.
func CodingSequenceCutoff(rbs *RibosomeBindingSite) interface{} {
	codingSequenceCutoff := 35

	lenProteinCodingSequence := len(rbs.ProteinCodingSequence)
	if codingSequenceCutoff > lenProteinCodingSequence {
		codingSequenceCutoff = lenProteinCodingSequence
	}
	return codingSequenceCutoff
}

// updateCodingSequence cuts off the coding sequence from the 5' to 3'
// direction to the value of the `codingSequenceCutoff` property value, and
// updates the `ProteinCodingSequence` field of `rbs`. NOTE: as principle,
// properties shouldn't update the fields of the `rbs` struct, but this is an
// initialization step so is an exception.
func updateCodingSequence(rbs *RibosomeBindingSite) interface{} {
	codingSequenceCutoff := rbs.PropertyValue(CodingSequenceCutoff).(int)
	rbs.ProteinCodingSequence = rbs.ProteinCodingSequence[:codingSequenceCutoff]
	return nil
}

/*****************************************************
Step 2: Calculate dG_start
*****************************************************/

// values taken from: https://pubs.acs.org/doi/10.1021/acssynbio.0c00394
var startCodonFreeEnergies map[string]float64 = map[string]float64{
	"AUG": -2.76,
	"GUG": -0.42,
	"UUG": 1.81,
	"CUG": 7.09,
	"AUC": 7.18,
	"AUA": 8.19,
	"AUU": 12.23,
	"CAU": 13.41,
	"GGA": 16.89,
	"UGC": 17.99,
	"CGC": 18.17,
	"UAG": 19.14,
}

// CDSStartCodonFreeEnergy (dG_start) is the free energy released when the
// start codon pairs to tRNA-fMet
func CDSStartCodonFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	cds := rbs.ProteinCodingSequence
	if len(cds) < 3 {
		return 0.0
	}
	startCodon := cds[:3]
	if dG_start, ok := startCodonFreeEnergies[startCodon]; ok {
		return dG_start
	} else {
		return 0.0
	}
}

/*****************************************************
Step 3: Find the Shine-Dalgarno Binding Site
*****************************************************/

const lenShineDalgarnoBindingSite int = 10

// shineDalgarnoBindingSite returns the shine-dalgarno binding site (of type
// `ShineDalgarnoBindingSite`) where the 16s rRNA of the organism binds to the
// binds to the mRNA sequence
func shineDalgarnoBindingSite(rbs *RibosomeBindingSite) interface{} {
	var bindingSites []ShineDalgarnoBindingSite

	lenFivePrimeUTR := len(rbs.FivePrimeUTR)
	// SD sequences can end at max 17 nucleotides away from the start codon
	maxAlignedDistFromStartCodonForSDSequece := 17
	sdSequenceSearchStartIdx := lenFivePrimeUTR - maxAlignedDistFromStartCodonForSDSequece
	if sdSequenceSearchStartIdx < lenShineDalgarnoBindingSite {
		sdSequenceSearchStartIdx = lenShineDalgarnoBindingSite
	}

	for bindingSiteEndPos := sdSequenceSearchStartIdx; bindingSiteEndPos <= lenFivePrimeUTR; bindingSiteEndPos++ {
		bindingSiteStartIdx := bindingSiteEndPos - lenShineDalgarnoBindingSite
		bindingSiteSDSequence := rbs.FivePrimeUTR[bindingSiteStartIdx:bindingSiteEndPos]
		sdBindingSite, ok := sdBindingSiteLookupTable[rbs.RibosomalRNA][bindingSiteSDSequence][rbs.Temperature]
		if ok && sdBindingSite.MRNAFivePrimeIdx != -1 {
			sdBindingSite.MessengerRNAStructure = sdBindingSite.MessengerRNAStructure[:sdBindingSite.MRNAThreePrimeIdx+1]
			sdBindingSite.MRNAFivePrimeIdx += bindingSiteStartIdx
			sdBindingSite.MRNAThreePrimeIdx += bindingSiteStartIdx

			nbNucleotidesBetweenMRNABindingSiteAndCDS := (lenFivePrimeUTR - 1) - sdBindingSite.MRNAThreePrimeIdx
			nbNucleotidesToCompleteASDSequence := sdBindingSite.RRNAFivePrimeIdx

			// the number of codons between the end of the aligned binding site and the
			// protein coding sequence
			sdBindingSite.AlignedSpacing = nbNucleotidesBetweenMRNABindingSiteAndCDS - nbNucleotidesToCompleteASDSequence

			// skip binding site where aligned spacing is 1 or less
			if sdBindingSite.AlignedSpacing <= 1 {
				continue
			}

			// compute upstream and downstream region energies
			upstreamRegion := rbs.FivePrimeUTR[:sdBindingSite.MRNAFivePrimeIdx]
			if len(upstreamRegion) != 0 {
				_, utrFreeEnergy := linearfold.ViennaRNAFold(upstreamRegion, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds, linearfold.DefaultBeamSize)
				sdBindingSite.UTRFreeEnergy = utrFreeEnergy
			}

			downstreamRegion := rbs.FivePrimeUTR[sdBindingSite.MRNAThreePrimeIdx+1:]
			if len(downstreamRegion) != 0 {
				_, dtrFreeEnergy := linearfold.ViennaRNAFold(downstreamRegion, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds, linearfold.DefaultBeamSize)
				sdBindingSite.DTRFreeEnergy = dtrFreeEnergy
			}

			// compute the free energy of the spacing region of this binding site
			sdBindingSite.AlignedSpacingFreeEnergy = compute_dG_spacing(sdBindingSite.AlignedSpacing)

			bindingSites = append(bindingSites, sdBindingSite)
		}
	}

	if len(bindingSites) > 0 {
		sort.Slice(bindingSites, func(i, j int) bool {
			iFreeEnergy := bindingSites[i].UTRFreeEnergy + bindingSites[i].SDFreeEnergy + bindingSites[i].DTRFreeEnergy + bindingSites[i].AlignedSpacingFreeEnergy
			jFreeEnergy := bindingSites[j].UTRFreeEnergy + bindingSites[j].SDFreeEnergy + bindingSites[j].DTRFreeEnergy + bindingSites[j].AlignedSpacingFreeEnergy
			if iFreeEnergy == jFreeEnergy {
				return bindingSites[i].MRNAThreePrimeIdx > bindingSites[j].MRNAThreePrimeIdx
			}
			return iFreeEnergy < jFreeEnergy
		})
		return &bindingSites[0]
	} else {
		return nil
	}
}

func HasShineDalgarnoBindingSite(rbs *RibosomeBindingSite) interface{} {
	switch rbs.PropertyValue(shineDalgarnoBindingSite).(type) {
	case nil:
		return false
	default:
		return true
	}
}

// Not used for further computation. Only for output.
func SDBindingSiteMRNAStructure(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *(rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite))
		return sdBindingSite.MessengerRNAStructure
	}
	return ""
}

// Not used for further computation. Only for output.
func SDBindingSiteRRNAStructure(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *(rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite))
		return sdBindingSite.RibosomalRNAStructure
	}
	return ""
}

// Not used for further computation. Only for output.
func BindingSiteMRNAStructure(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *(rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite))
		return sdBindingSite.MessengerRNA
	}
	return ""
}

// Not used for further computation. Only for output.
func LenSDSequence(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *(rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite))
		return sdBindingSite.MRNAThreePrimeIdx - sdBindingSite.MRNAFivePrimeIdx + 1
	} else {
		return 0
	}
}

func BindingSiteFivePrimeIdx(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite)
		return sdBindingSite.MRNAFivePrimeIdx
	} else {
		return -1
	}
}

// Not used for further computation. Only for output.
func BindingSiteThreePrimeIdx(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite)
		return sdBindingSite.MRNAThreePrimeIdx
	} else {
		return -1
	}
}

/*****************************************************
Step 4: Calculate dG_spacing
*****************************************************/

// SpacingRegionFreeEnergy (dG_spacing) is the free energy penalty for the
// stretching / compression of the spacing region between the shine-dalgarno
// binding site and start codon.
func SpacingRegionFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite)
		return sdBindingSite.AlignedSpacingFreeEnergy
	} else {
		return 0.0
	}
}

// Not used for further computation. Only for output.
func AlignedSpacing(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite)
		return sdBindingSite.AlignedSpacing
	} else {
		return 0
	}
}

// sb0c00394_si_001.pdf page 13
func compute_dG_spacing(alignedSpacing int) float64 {
	// var optimalSpacing int
	// if organism == "Bacillus subtilis subsp. subtilis str. 168" {
	// 	optimalSpacing = 6
	// } else {
	// }
	optimalSpacing := 5

	alignedSpacingOptDiff := float64(alignedSpacing - optimalSpacing)
	if alignedSpacingOptDiff > 0.0 {
		// RBS Calculator v1 - v2.0 parameters:
		c1 := 0.048 // kcal/mol/nt^2
		c2 := 0.24  // kcal/mol/nt

		// RBS Calculator v2.1 parameters:
		// c1 := 0.136 // kcal/mol/nt2
		// c2 := 0.068 // kcal/mol/nt
		return c1*alignedSpacingOptDiff*alignedSpacingOptDiff + c2*alignedSpacingOptDiff
	} else if alignedSpacingOptDiff < 0.0 {
		// RBS Calculator v1 - v2.0 parameters:
		c3 := 12.2 // kcal/mol
		c4 := 2.5  // nt^(-1)
		c5 := 2.0  // nt
		c6 := 3.0

		// RBS Calculator v2.1 parameters:
		// c3 := 12.9  // kcal/mol
		// c4 := 1.15  // nt-1
		// c5 := 0.154 // nt
		// c6 := 4.5
		return c3 / math.Pow(1+math.Exp(c4*(alignedSpacingOptDiff+c5)), c6)
	} else {
		return 0.0
	}
}

/*****************************************************
Step 5: Calculate the free energy of the mRNA sequence (dG_mRNA)
*****************************************************/

// UsedMRNASequence is the mRNA sequence used to calculate the dG_mRNA term of
// the model
func UsedMRNASequence(rbs *RibosomeBindingSite) interface{} {
	return rbs.FivePrimeUTR + rbs.ProteinCodingSequence
}

// UsedMRNADotBracketStructure is the folded structure of the used mRNA sequence.
func UsedMRNADotBracketStructure(rbs *RibosomeBindingSite) interface{} {
	mRNA := rbs.PropertyValue(UsedMRNASequence).(string)
	structure, _ := linearfold.ViennaRNAFold(mRNA, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds, linearfold.DefaultBeamSize)
	return structure
}

// MRNAFreeEnergy (dG_mRNA) is the free energy of the folded mRNA sequence.
func MRNAFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	mRNA := rbs.PropertyValue(UsedMRNASequence).(string)
	mRNAStructure := rbs.PropertyValue(UsedMRNADotBracketStructure).(string)

	dG_mRNA, _, err := mfe.MinimumFreeEnergy(mRNA, mRNAStructure, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds)
	if err != nil {
		panic(err)
	}
	return dG_mRNA
}

/*****************************************************
Step 6: Calculate dG_mRNA_rRNA

dG_mRNA_rRNA is the hybridization energy of the between the mRNA and 16s rRNA
of the 30S ribosomal subunit of the host organism together with other mRNA-mRNA
interactions do not require unfolding during translation initiation.

To calculate this term, we have to:
1. Calculate the free energy of the mRNA-mRNA interactions that occur before the
	 shine-dalgarno binding site (the site where the 16s rRNA binds to the mRNA) -
	 dG_pre_ribosome.
2. Calculate the free energy of the hybridization between the 16s rRNA and the
	 mRNA at the shine-dalgarno binding site - dG_SD_aSD.
3. Calculate the free energy of the mRNA-mRNA interactions that occur after the
	 shine-dalgarno binding site - dG_post_ribosome.

dG_mRNA_rRNA is the sum of these three terms.
*****************************************************/

/**********************************
Calculate dG_pre_ribosome
**********************************/

// MRNAPreRibosomeSequence is the sequence of the 5' UTR till where the
// ribosome binds to the mRNA
func MRNAPreRibosomeSequence(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		bindingSiteFivePrime := rbs.PropertyValue(BindingSiteFivePrimeIdx).(int)
		return rbs.FivePrimeUTR[:bindingSiteFivePrime]
	} else {
		// rbs doesn't have sd binding site so return the whole five prime UTR
		return rbs.FivePrimeUTR
	}
}

// MRNAPreRibosomeDotBracketStructure is the folded dot-bracket structure
// of the `mrnaPreSDBindingSite` sequence
func MRNAPreRibosomeDotBracketStructure(rbs *RibosomeBindingSite) interface{} {
	preRibosomeMRNA := rbs.PropertyValue(MRNAPreRibosomeSequence).(string)
	if preRibosomeMRNA == "" {
		return ""
	} else {
		preRibosomeMRNAStructure, _ := linearfold.ViennaRNAFold(preRibosomeMRNA, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds, linearfold.DefaultBeamSize)
		return preRibosomeMRNAStructure
	}
}

// PreRibosomeFreeEnergy (dG_pre_ribosome) is the free energy of the mRNA-mRNA
// interactions in the sequence before the shine-dalgarno binding site
func PreRibosomeFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	mrnaPreRibosomeDotBracketStructure := rbs.PropertyValue(MRNAPreRibosomeDotBracketStructure).(string)
	if mrnaPreRibosomeDotBracketStructure == "" {
		return 0.0
	} else {
		mrnaPreRibosomeSequence := rbs.PropertyValue(MRNAPreRibosomeSequence).(string)
		dG_pre_ribosome, _, err := mfe.MinimumFreeEnergy(mrnaPreRibosomeSequence, mrnaPreRibosomeDotBracketStructure, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds)
		if err != nil {
			panic(err)
		}
		return dG_pre_ribosome
	}
}

/**********************************
Calculate dG_SD_aSD
**********************************/

// ShineDalgarnoHybridizationFreeEnergy (dG_SD_aSD) is the hybridization free
// energy that occurs at the shine-dalgarno binding site (the site where the
// 16s rRNA of the 30S ribosomal subunit of the host organism binds to the mRNA
// strand)
func ShineDalgarnoHybridizationFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := *rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite)
		// add the duplex free energy adjustment
		dG_SD_aSD := sdBindingSite.SDFreeEnergy + duplexFreeEnergyAdjustment(rbs.Temperature)
		return dG_SD_aSD
	} else {
		return 0.0
	}
}

// Taken from:
// https://github.com/reisalex/SynBioMTS/blob/d7d46af7981f90d2147cbc0463c5254b16c5c92c/examples/RBS/models/PyVRNA.py#L248
// duplexFreeEnergyAdjustment is the free energy adjustment added when we
// compute the hybridization energy between two RNA strands.
func duplexFreeEnergyAdjustment(temperature float64) float64 {
	kB := 0.00198717 // Boltzmann constant in kcal/mol/K
	a := []float64{-3.983035, 301.797, 522528.9, 69.34881, 999.974950}
	pH2O := a[4] * (1 - math.Pow(temperature+a[0], 2.0)*(temperature+a[1])/a[2]/(temperature+a[3])) / 18.0152

	return -kB * (temperature + 273.15) * math.Log(pH2O)
}

/**********************************
Calculate dG_post_ribosome
**********************************/

// MRNAPostRibosomeSequence is the mRNA sequence after the footprint of the
// ribosome
func MRNAPostRibosomeSequence(rbs *RibosomeBindingSite) interface{} {
	cds := rbs.ProteinCodingSequence
	// as per https://pubmed.ncbi.nlm.nih.gov/28158713/, length of the ribosome
	// footprint is 13 nucleotides from the start codon, so we reduce start the
	// coding sequence from its 13th nucleotide
	postRibosomeStartIdx := 13
	lenCDS := len(cds)
	if lenCDS < postRibosomeStartIdx+1 {
		postRibosomeStartIdx = lenCDS
	}
	return cds[postRibosomeStartIdx:]
}

// MRNAPostRibosomeDotBracketStructure is the folded dot-bracket structure of
// `mrnaPostRibosomeSequence`
func MRNAPostRibosomeDotBracketStructure(rbs *RibosomeBindingSite) interface{} {
	mrnaPostRibosomeSequence := rbs.PropertyValue(MRNAPostRibosomeSequence).(string)
	if mrnaPostRibosomeSequence == "" {
		return ""
	} else {
		dotBracketStructure, _ := linearfold.ViennaRNAFold(mrnaPostRibosomeSequence, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds, linearfold.DefaultBeamSize)
		return dotBracketStructure
	}
}

// PostRibosomeFreeEnergy (dG_post_ribosome) is the free energy of the mRNA-mRNA
// interactions in the sequence after the ribosome
func PostRibosomeFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	mrnaPostRibosomeDotBracketStructure := rbs.PropertyValue(MRNAPostRibosomeDotBracketStructure).(string)
	if mrnaPostRibosomeDotBracketStructure == "" {
		return 0.0
	} else {
		mrnaPostRibosomeSequence := rbs.PropertyValue(MRNAPostRibosomeSequence).(string)
		dG_pre_ribosome, _, err := mfe.MinimumFreeEnergy(mrnaPostRibosomeSequence, mrnaPostRibosomeDotBracketStructure, rbs.Temperature, energy_params.Andronescu2007, mfe.NoDanglingEnds)
		if err != nil {
			panic(err)
		}
		return dG_pre_ribosome
	}
}

/**********************************
Finally, compute dG_mRNA_rRNA
**********************************/

// MRNARRNAHybridizationFreeEnergy is the dG_mRNA_rRNA term of the model
func MRNARRNAHybridizationFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	dG_pre_ribosome := rbs.PropertyValue(PreRibosomeFreeEnergy).(float64)
	dG_post_ribosome := rbs.PropertyValue(PostRibosomeFreeEnergy).(float64)
	dG_SD_aSD := rbs.PropertyValue(ShineDalgarnoHybridizationFreeEnergy).(float64)
	return dG_pre_ribosome + dG_post_ribosome + dG_SD_aSD
}

/*****************************************************
Step 7: Calculate dG_standby
*****************************************************/

// MRNAPreRibosomeSecondaryStructure is the `SecondaryStructure` of the
// pre-ribosome mRNA sequence
func MRNAPreRibosomeSecondaryStructure(rbs *RibosomeBindingSite) interface{} {
	mrnaPreRibosomeDotBracketStructure := rbs.PropertyValue(MRNAPreRibosomeDotBracketStructure).(string)
	_, secondaryStructure, err := SecondaryStructureFromDotBracket(mrnaPreRibosomeDotBracketStructure)
	if err != nil {
		panic(err)
	}
	return secondaryStructure
}

// ribosomes likely bind to upstream standby sites. However, upstream standby
// sites may have structures that cause them to be less favorable and hence
// have a positive free energy associated with the site. HairpinModule returns
// the most favorable upstream standby site.
func HairpinModule(rbs *RibosomeBindingSite) interface{} {
	mrnaPreRibosomeSecondaryStructure := rbs.PropertyValue(MRNAPreRibosomeSecondaryStructure).(*SecondaryStructure)
	if mrnaPreRibosomeSecondaryStructure == nil {
		return nil
	}

	// get all the upstream hairpin modules
	hairpinModules := hm.HairpinModules(mrnaPreRibosomeSecondaryStructure.Structures)

	// for each hairpin module, add the sliding free energy penalty (the amount of
	// energy required to unfold the downstream hairpin modules
	hairpinModules = hm.AddDownstreamModulesSlidingFreeEnergy(hairpinModules)

	if len(hairpinModules) > 0 {
		// choose the hairpin module with the least standby free energy. If
		// the standby energies are the same, choose the hairpin module closer to
		// the shine-dalgarno binding site
		sort.Slice(hairpinModules, func(i, j int) bool {
			if hairpinModules[i].StandbyFreeEnergy == hairpinModules[j].StandbyFreeEnergy {
				return hairpinModules[i].Hairpin.SingleStrandedThreePrimeIdx > hairpinModules[j].Hairpin.SingleStrandedThreePrimeIdx
			}
			return hairpinModules[i].StandbyFreeEnergy < hairpinModules[j].StandbyFreeEnergy
		})
		return hairpinModules[0]
	} else {
		return nil
	}
}

func HasHairpinModule(rbs *RibosomeBindingSite) interface{} {
	switch rbs.PropertyValue(HairpinModule).(type) {
	case nil:
		return false
	default:
		return true
	}
}

// StandbyModuleTotalFreeEnergy (dG_standby) is the sum of dG_unfolding, dG_sliding
// and dG_distortion
func StandbyModuleTotalFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.StandbyFreeEnergy
	} else {
		return 0.0
	}
}

// StandbyModuleUnfoldingFreeEnergy is the dG_unfolding term of the model.
// Not used for further computation. Only for output.
func StandbyModuleUnfoldingFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.UnfoldingFreeEnergy
	} else {
		return 0.0
	}
}

// StandbyModuleSlidingFreeEnergy is the dG_sliding term of the model.
// Not used for further computation. Only for output.
func StandbyModuleSlidingFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.SlidingFreeEnergy
	} else {
		return 0.0
	}
}

// StandbyModuleDistortionFreeEnergy is the dG_distortion term of the model.
// Not used for further computation. Only for output.
func StandbyModuleDistortionFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.DistortionFreeEnergy
	} else {
		return 0.0
	}
}

// Not used for further computation. Only for output.
func HairpinModuleDistalLength(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.CalculateDistalLength()
	} else {
		return 0
	}
}

// Not used for further computation. Only for output.
func HairpinModuleProximalLength(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.CalculateProximalLength()
	} else {
		return 0
	}
}

// Not used for further computation. Only for output.
func HairpinModuleHairpinHeight(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.HairpinHeight
	} else {
		return 0
	}
}

// Not used for further computation. Only for output.
func HairpinModuleSurfaceArea(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.CalculateSurfaceArea()
	} else {
		return 0
	}
}

// Not used for further computation. Only for output.
func HairpinModuleFivePrimeIdx(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.FivePrimeIdx()
	} else {
		return -1
	}
}

// Not used for further computation. Only for output.
func HairpinModuleThreePrimeIdx(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasHairpinModule).(bool) {
		hairpinModule := rbs.PropertyValue(HairpinModule).(*hm.HairpinModule)
		return hairpinModule.ThreePrimeIdx()
	} else {
		return -1
	}
}

/*****************************************************
Step 8: Calculate dG_stack
*****************************************************/

// SpacingSequence is the nucleotide sequence which separates the end of the
// shine-dalgarno binding site from the start of the coding sequence
func SpacingSequence(rbs *RibosomeBindingSite) interface{} {
	if rbs.PropertyValue(HasShineDalgarnoBindingSite).(bool) {
		sdBindingSite := rbs.PropertyValue(shineDalgarnoBindingSite).(*ShineDalgarnoBindingSite)
		spacingRegionFivePrimeIdx := sdBindingSite.MRNAThreePrimeIdx + 1
		spacingRegion := rbs.FivePrimeUTR[spacingRegionFivePrimeIdx:]
		return spacingRegion
	} else {
		return ""
	}
}

// LenSpacingSequence is the length of `spacingSequence`
func LenSpacingSequence(rbs *RibosomeBindingSite) interface{} {
	return len(rbs.PropertyValue(SpacingSequence).(string))
}

// stackingEnergy is the
var stackingEnergy = map[byte]map[byte]float64{
	'A': map[byte]float64{
		'A': 0.0451,
		'U': 0.1282,
		'C': -0.0558,
		'G': 0.0451,
	},
	'U': map[byte]float64{
		'A': 0.1282,
		'U': 0.2603,
		'C': 0.0518,
		'G': 0.1282,
	},
	'C': map[byte]float64{
		'A': -0.0558,
		'U': 0.0518,
		'C': -0.1568,
		'G': -0.0558,
	},
	'G': map[byte]float64{
		'A': 0.0451,
		'U': 0.1282,
		'C': -0.0558,
		'G': 0.0451,
	},
}

// SpacingSequenceStackFreeEnergy is the dG_stack term of the model
func SpacingSequenceStackFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	spacingSequence := rbs.PropertyValue(SpacingSequence).(string)
	dG_stack := 0.0
	for i := 1; i < len(spacingSequence); i++ {
		currNucleotide := spacingSequence[i]
		prevNucleotide := spacingSequence[i-1]
		dG_stack += stackingEnergy[currNucleotide][prevNucleotide]
	}
	return dG_stack
}

// TotalFreeEnergy (dG_total) is the total free energy of the ribosome binding site
func TotalFreeEnergy(rbs *RibosomeBindingSite) interface{} {
	dG_spacing := rbs.PropertyValue(SpacingRegionFreeEnergy).(float64)
	dG_start := rbs.PropertyValue(CDSStartCodonFreeEnergy).(float64)
	dG_standby := rbs.PropertyValue(StandbyModuleTotalFreeEnergy).(float64)
	dG_mRNA_rRNA := rbs.PropertyValue(MRNARRNAHybridizationFreeEnergy).(float64)
	dG_mRNA := rbs.PropertyValue(MRNAFreeEnergy).(float64)
	dG_stack := rbs.PropertyValue(SpacingSequenceStackFreeEnergy).(float64)
	dG_total := dG_standby + dG_mRNA_rRNA + dG_spacing + dG_start + dG_stack - dG_mRNA
	return dG_total
}

/********************
Helper functions used to compute properties
*********************/

func TranslationInitiationRate(rbs RibosomeBindingSite) (translationInitiationRate float64) {
	totalFreeEnergy := rbs.PropertyValue(TotalFreeEnergy).(float64)
	// slope, intercept := -0.053698259, 9.329199575
	// slope, intercept := -0.05490834402441548, 9.353335814085765
	// slope, intercept := -0.05369825930168624, 9.329199575270719
	slope, intercept := -0.45, 7.11720550316
	// slope, intercept := -0.0536982593017, 9.32919957527
	translationInitiationRate = math.Exp(totalFreeEnergy*slope + intercept)
	// translationInitiationRate = totalFreeEnergy*slope + intercept
	return
}
