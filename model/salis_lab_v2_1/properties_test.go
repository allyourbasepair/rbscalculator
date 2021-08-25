// +build gen_csv

// Since we don't want the tests below to run when `go test ./...` is used, the
// above build contrain flag ensures the tests below are only run when we pass
// the `gen_csv` tag to go test with the tags argument
// (`go test ./... -tags gen_csv`).
package salis_lab_v2_1

import (
	"github.com/allyourbasepair/rbscalculator/model"
)

// `cd` into this directory and run with
// `go test -timeout 0 -run ^ExampleComputeProperties_forTrainDataset`
func ExampleComputeProperties_forTrainDataset() {
	datasetName := "train"

	idColIdx, datasetColIdx, organismColIdx, tempColIdx, proteinColIdx, fivePrimeUTRColIdx, cdsColNum, proteinMeanColIdx, proteinStdColIdx, ribosomalRNAColIdx := 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

	var otherInformationColIdxNameMap map[int]string = map[int]string{
		idColIdx:          "ID",
		datasetColIdx:     "Dataset",
		organismColIdx:    "Organism",
		proteinColIdx:     "Protein",
		proteinMeanColIdx: "Protein Mean",
		proteinStdColIdx:  "Protein Std",
	}

	propertiesToOutputToCSV := []model.RBSPropertyFunc{
		// Step 1: Initialization of `rbs` struct
		FivePrimeUTRCutoff,
		CodingSequenceCutoff,

		// Step 2: Calculate dG_start
		CDSStartCodonFreeEnergy,

		// Step 3: Find the Shine-Dalgarno Binding Site
		HasShineDalgarnoBindingSite,
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
		MRNAPreRibosomeDotBracketStructure,
		PreRibosomeFreeEnergy,

		// Calculate dG_SD_aSD
		ShineDalgarnoHybridizationFreeEnergy,

		// Calculate dG_post_ribosome
		MRNAPostRibosomeDotBracketStructure,
		PostRibosomeFreeEnergy,

		MRNARRNAHybridizationFreeEnergy,

		// Step 7: Calculate dG_standby
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
		SpacingSequenceStackingFreeEnergy,

		TotalFreeEnergy,
	}

	model.ComputePropertiesForDataset(datasetName, PropertiesToCompute, propertiesToOutputToCSV, fivePrimeUTRColIdx, cdsColNum, tempColIdx, ribosomalRNAColIdx, otherInformationColIdxNameMap)
}
