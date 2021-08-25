package rbs_calculator

import (
	"fmt"
	"os"
	"regexp"
	"sort"

	"github.com/allyourbasepair/rbscalculator/model"
	rbs_model "github.com/allyourbasepair/rbscalculator/model/salis_lab_v2_1"
	"github.com/olekukonko/tablewriter"
)

/*****************************************************
June, 5, 2021
The Ribosome Binding Site Calculator

Translation is the process by which a messenger RNA (mRNA) is 'decoded' to a
protein by a ribosome.

When trying to synthesize proteins of our choice, we have to design a mRNA
sequence with a  5' untranslated region, the protein coding sequence and a
terminator sequence. The ribosome's interaction with the mRNA sequence
determines the amount of the protein that will be synthesized. Thus, we need
a model to understand these interactions and how they relate to the amount of
protein produced.

The ribosome binding site calculator returns a list of potential places (binding
sites) a ribosome might bind to a mRNA. Each binding site breaks down a mRNA
sequence into a 5' untranslated region and protein coding sequence, and includes
information about the predicted translation rate of that binding site. This is
useful because it allows one to examine where a ribosome is likely to bind to
the mRNA sequence by examining the relative differences between the translation
initiation rates of each binding site.

The other intended use of this calculator is to help one increase synthesis of
a desired protein by the user re-desinging their mRNA sequence to increase
translation initiation rate which *should* lead to an increase in the amount of
the protein produced in-vitro.

How the calculator works:
To calculate the translation initiation rate of a binding site, we need a model
of a ribosome's interactions with a mRNA strand and the relationship of the
interactions with the translation initiation rate of the coding sequence of the
mRNA strand.

The current model of this RBS calculator is based upon the research done by
The Salis Lab at The Penn State University. For more information on this model,
please have a look at the `properties.go` file in the `salis_lab_v2_1`
subpackage (./model/salis_lab_v2_1/properties.go).

For information of how to develop a RBS Calculator of your own, please
have a look at the `model.go` file in the `model` subpackage (./model/model.go).

*****************************************************/

// StartCodon specifies the start codon of the protein coding sequence of a
// mRNA strand.
type StartCodon string

const (
	AUG StartCodon = "AUG"
	GUG            = "GUG"
	UUG            = "UUG"
	CUG            = "CUG"
	AUC            = "AUC"
	AUA            = "AUA"
	AUU            = "AUU"
	CAU            = "CAU"
	GGA            = "GGA"
	UGC            = "UGC"
	CGC            = "CGC"
	UAG            = "UAG"
)

// Organism16SrRNAMap is a map of an organism to its 16S ribosomal RNA.
var Organism16SrRNAMap map[string]string = map[string]string{
	"Bacillus subtilis subsp. subtilis str. 168":                       "ACCUCCUUU",
	"Bacteroides thetaiotaomicron VPI-5482":                            "ACCUCCUUA",
	"Corynebacterium glutamicum B-2784":                                "ACCUCCUUU",
	"Escherichia coli BL21(DE3)":                                       "ACCUCCUUA",
	"Escherichia coli str. K-12 substr. DH10B":                         "ACCUCCUUA",
	"Escherichia coli str. K-12 substr. MG1655":                        "ACCUCCUUA",
	"Pseudomonas fluorescens A506":                                     "ACCUCCUUA",
	"Salmonella enterica subsp. enterica serovar Typhimurium str. LT2": "ACCUCCUUA",
}

// DefaultStartCodons are AUG, GUG and CUG
var DefaultStartCodons []StartCodon = []StartCodon{AUG, GUG, CUG}

// RibosomeBindingSites returns a list of ribosome binding sites for a given
// messenger RNA (mRNA) strand and 16s rRNA. Each site contains information
// about the five prime untranslated region, the protein coding sequence, the
// translation initiation rate as well as all the properties that were computed
// to figure out the translation initiation rate.
//
// The translation initiation rates of the binding sites can be compared with
// each other to help you figure out where a ribosome is likely to bind to
// your mRNA strand.
//
// The output can also be used to increase the amount of synthesis of your
// desired protein. Redesign your mRNA sequence to have higher translation
// initiation rates and this *should* lead to an increase in the amount
// of protein synthesized by your mRNA strand in-vitro.
//
// Use the exported map `Organism16SrRNAMap` to find the 16s rRNA of the
// organism the mRNA strand is inserted into.
func RibosomeBindingSites(ribosomalRNA, mRNA string, temperatureInCelsius float64, startCodons []StartCodon) (ribosomeBindingSites []model.RibosomeBindingSite) {
	ribosomalRNA = model.CleanRNA(ribosomalRNA)
	mRNA = model.CleanRNA(mRNA)

	// for each start codon, find all occurrences of it in the mRNA sequence,
	// create a `model.RibosomeBindingSite` struct for it and compute the
	// translation initiation rate for each struct
	for _, startCondon := range startCodons {
		startCodonRegex := regexp.MustCompile(string(startCondon))
		matches := startCodonRegex.FindAllStringIndex(mRNA, -1)
		for _, match := range matches {
			startCodonIdx := match[0]
			fivePrimeUTR, codingSequence := mRNA[:startCodonIdx], mRNA[startCodonIdx:]

			_, ribosomeBindingSite := TranslationInitiationRate(fivePrimeUTR, codingSequence, ribosomalRNA, temperatureInCelsius)
			ribosomeBindingSites = append(ribosomeBindingSites, ribosomeBindingSite)
		}
	}

	// sort by length of 5' UTR (start position) in ascending order
	sort.Slice(ribosomeBindingSites, func(i, j int) bool {
		return ribosomeBindingSites[i].PropertyValue(model.StartPosition).(int) < ribosomeBindingSites[j].PropertyValue(model.StartPosition).(int)
	})
	return
}

// SortByTranslationInitiationRate sorts a list of binding sites by their
// translation initiation rates in descending order
func SortByTranslationInitiationRate(ribosomeBindingSites []model.RibosomeBindingSite) []model.RibosomeBindingSite {
	sort.Slice(ribosomeBindingSites, func(i, j int) bool {
		if ribosomeBindingSites[i].TranslationInitiationRate == ribosomeBindingSites[j].TranslationInitiationRate {
			return len(ribosomeBindingSites[i].FivePrimeUTR) < len(ribosomeBindingSites[j].FivePrimeUTR)
		}
		return ribosomeBindingSites[i].TranslationInitiationRate > ribosomeBindingSites[j].TranslationInitiationRate
	})
	return ribosomeBindingSites
}

// TranslationInitiationRate returns the the translation initiation rate of a
// ribsome binding site as well as the binding site (as a
// `model.RibosomeBindingSite` struct) with the properties computed to calculate
// the translation initiation rate
func TranslationInitiationRate(fivePrimeUTR, proteinCodingSequence, ribosomalRNA string, temperateureInCelsius float64) (translationInitiationRate float64, bindingSiteWithProperties model.RibosomeBindingSite) {
	rbs := model.RibosomeBindingSite{
		FivePrimeUTR:          fivePrimeUTR,
		ProteinCodingSequence: proteinCodingSequence,
		Temperature:           temperateureInCelsius,
		RibosomalRNA:          ribosomalRNA,
		Properties:            make(map[string]interface{}),
		// we use inf as a sanity check to ensure `TranslationInitiationRate` has
		// been computed for the binding site
		TranslationInitiationRate: inf,
	}

	// compute the properties of the interactions between the mRNA strand
	// and ribosome required to calculate the translation initiation rate
	rbs.ComputeProperties(model.DefaultPropertiesToComputeBefore)
	rbs.ComputeProperties(rbs_model.PropertiesToCompute)

	// calculate the translation initiation rate
	translationInitiationRate = rbs_model.ComputeTranslationInitiationRate(rbs)
	rbs.TranslationInitiationRate = translationInitiationRate

	rbs.ComputeProperties(model.DefaultPropertiesToComputeAfter)

	if rbs.TranslationInitiationRate == inf {
		panic("failed to calculate the translation initiation rate for the mRNA sequence. ensure the `TranslationInitiationRate` field of the `model.RibosomeBindingSite` struct is set by the `rbs_model.TranslationInitiationRate` func.")
	}

	return translationInitiationRate, rbs
}

// PrintBindingSites prints the important properties of a binding site.
// The optional argument `additionalPropertiesToPrint` specifies the computed
// properties of the binding site that will be printed.
func PrintBindingSites(bindingSites []model.RibosomeBindingSite, includeSequences, includeStructures bool, additionalPropertiesToPrint ...PropertyToPrint) {
	table := tablewriter.NewWriter(os.Stdout)
	table.SetRowLine(true)
	table.SetAutoFormatHeaders(false)

	// add the basic important properties of the binding site
	propertiesToPrint := []PropertyToPrint{
		{property: model.StartPosition, columnHeader: "Start position"},
		{property: model.TranslationInitiationRate, columnHeader: "TIR"},
	}

	// add the additional properties
	propertiesToPrint = append(propertiesToPrint, additionalPropertiesToPrint...)

	if includeSequences {
		sequenceProperties := []PropertyToPrint{
			{property: model.FivePrimeUTR, columnHeader: "5' Untranslated Region"},
			{property: model.ProteinCodingSequence, columnHeader: "Protein Coding Sequence"},
		}

		propertiesToPrint = append(propertiesToPrint, sequenceProperties...)
	}

	if includeStructures {
		propertiesToPrint = append(propertiesToPrint, salisLabRBSModelStructureProperties...)
	}

	table.SetHeader(propertiesToPrintHeader(propertiesToPrint))

	// add the relevant information for each binding site
	for _, bindingSite := range bindingSites {
		bindingSiteInfo := rbsPropertyValues(bindingSite, propertiesToPrint)
		table.Append(bindingSiteInfo)
	}

	// finally, render the table
	table.Render()
}

// propertiesToPrintHeader returns a list of the `columnHeader` fields of each
// `PropertyToPrint` struct
func propertiesToPrintHeader(properties []PropertyToPrint) (ret []string) {
	for _, propertyToPrint := range properties {
		ret = append(ret, propertyToPrint.columnHeader)
	}
	return
}

// rbsPropertyValues returns a list of the values of the `property` fields
// of each `PropertyToPrint` struct
func rbsPropertyValues(rbs model.RibosomeBindingSite, properties []PropertyToPrint) (ret []string) {
	for _, propertyToPrint := range properties {
		ret = append(ret, toString(rbs.PropertyValue(propertyToPrint.property)))
	}
	return
}

// toString returns the string value of `i`
func toString(i interface{}) string {
	return fmt.Sprint(i)
}

// PropertyToPrint specifies the computed property of a ribosome binding site
// that will be printed when using the `PrintBindingSites` func
type PropertyToPrint struct {
	property     model.RBSPropertyFunc
	columnHeader string
}

var inf float64 = 1000000000.0

// salisLabRBSModelStructureProperties are all the computed properties in the
// Salis Lab v2.1 RBS Calc model that contain the dot-bracket structures
// of the initial and final structure states.
//
// This variable should ideally be in the `salis_lab_v2_1` subpackage. However,
// as the `PropertyToPrint` struct is defined in this package, adding this
// variable to the `salis_lab_v2_1` package would lead to a cyclical import
// which is not allowed in Go.
var salisLabRBSModelStructureProperties []PropertyToPrint = []PropertyToPrint{
	{property: rbs_model.UsedMRNADotBracketStructure, columnHeader: "Initial state"},
	{property: rbs_model.MRNAPreRibosomeDotBracketStructure, columnHeader: "Final state (pre ribosome)"},
	{property: rbs_model.SDBindingSiteMRNAStructure, columnHeader: "Final state (mRNA shine dalgarno binding site)"},
	{property: rbs_model.SpacingSequenceDotBracketStructure, columnHeader: "Final state (spacing)"},
	{property: rbs_model.RibosomeFootprintDotBracketStructure, columnHeader: "Final state (ribosome footprint)"},
	{property: rbs_model.MRNAPostRibosomeDotBracketStructure, columnHeader: "Final state (post ribosome)"},
	{property: rbs_model.SDBindingSiteRRNAStructure, columnHeader: "Final state (16S rRNA shine dalgarno binding site)"},
	{property: rbs_model.FinalStateDotBracketStructure, columnHeader: "Full final state"},
}
