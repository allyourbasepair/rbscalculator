package hairpin_modules

import (
	"math"

	. "github.com/TimothyStiles/poly/secondary_structure"
)

type HairpinModule struct {
	HairpinHeight                                int
	distalSingleStranded, proximalSingleStranded *SingleStrandedRegion
	Hairpin                                      Hairpin
	downstreamHairpinModules                     []*HairpinModule
	StandbyFreeEnergyWithoutSliding              float64
	StandbyFreeEnergy                            float64
	SlidingFreeEnergy                            float64
	DistortionFreeEnergy                         float64
	UnfoldingFreeEnergy                          float64
}

func (hairpinModule *HairpinModule) CalculateDistalLength() int {
	distalSingleStranded := hairpinModule.distalSingleStranded
	if distalSingleStranded != nil {
		length := distalSingleStranded.ThreePrimeIdx - distalSingleStranded.FivePrimeIdx + 1
		return length
	} else {
		return 0
	}
}

func (hairpinModule *HairpinModule) CalculateProximalLength() int {
	proximalSingleStranded := hairpinModule.proximalSingleStranded
	if proximalSingleStranded != nil {
		length := proximalSingleStranded.ThreePrimeIdx - proximalSingleStranded.FivePrimeIdx + 1
		return length
	} else {
		return 0
	}
}

const maxHairpinHeight = 15

// hairpins have a maximum height of 15
// height is calculated as the length of a Hairpin (divided by 2)
// Thus, it includes the height of one side of the stem + half of the nucleotides
// in the single stranded region
func (hairpinModule *HairpinModule) calculateHairpinHeight() int {
	hairpin := hairpinModule.Hairpin
	hairpinHeight := int(math.Ceil(float64(hairpin.Stem.ClosingThreePrimeIdx-hairpin.Stem.ClosingFivePrimeIdx+1) / 2.0))
	if hairpinHeight > maxHairpinHeight {
		hairpinHeight = maxHairpinHeight
	}
	return hairpinHeight
}

func newHairpinModule(hairpin *Hairpin, distalSingleStranded, proximalSingleStranded *SingleStrandedRegion) HairpinModule {
	hairpinModule := HairpinModule{
		distalSingleStranded:   distalSingleStranded,
		proximalSingleStranded: proximalSingleStranded,
		Hairpin:                *hairpin,
	}
	hairpinModule.optimizeHairpin()
	hairpinModule.HairpinHeight = hairpinModule.calculateHairpinHeight()
	return hairpinModule
}

func (hairpinModule *HairpinModule) CalculateSurfaceArea() int {
	distalLength := hairpinModule.CalculateDistalLength()
	hairpinHeight := hairpinModule.calculateHairpinHeight()
	proximalLength := hairpinModule.CalculateProximalLength()
	return maxHairpinHeight + distalLength + proximalLength - hairpinHeight
}

func calculateFreeEnergyDistortion(hairpinModule HairpinModule) float64 {
	return dG_distortion(hairpinModule.CalculateSurfaceArea())
}

func (hairpinModule *HairpinModule) optimizeHairpin() {
	// set best module to current module
	bestHairpinModule := *hairpinModule
	bestHairpinModule.DistortionFreeEnergy = calculateFreeEnergyDistortion(bestHairpinModule)
	bestFreeEnergyStandby := bestHairpinModule.DistortionFreeEnergy

	// get the next module
	nextHairpinModule, dG_unfolding, didUnfold := unfoldHairpin(bestHairpinModule)

	var currUnfoldingFreeEnergy float64
	// if we could unfold the hairpin in the module, continue
	for didUnfold {
		// update the current dG_unfolding
		currUnfoldingFreeEnergy += dG_unfolding
		// calculate the standby energy and add the unfolding term
		nextHairpinModule.DistortionFreeEnergy = calculateFreeEnergyDistortion(nextHairpinModule)
		nextHairpinModule.UnfoldingFreeEnergy = currUnfoldingFreeEnergy
		nextFreeEnergyStandby := nextHairpinModule.DistortionFreeEnergy + nextHairpinModule.UnfoldingFreeEnergy

		if nextFreeEnergyStandby < bestFreeEnergyStandby {
			bestHairpinModule = nextHairpinModule
			bestFreeEnergyStandby = nextFreeEnergyStandby
			nextHairpinModule, dG_unfolding, didUnfold = unfoldHairpin(bestHairpinModule)
		} else {
			break
		}
	}

	// no more to unfold in stem so return current best hairpin module
	hairpinModule.distalSingleStranded = bestHairpinModule.distalSingleStranded
	hairpinModule.proximalSingleStranded = bestHairpinModule.proximalSingleStranded
	hairpinModule.Hairpin = bestHairpinModule.Hairpin
	hairpinModule.DistortionFreeEnergy = bestHairpinModule.DistortionFreeEnergy
	hairpinModule.UnfoldingFreeEnergy = bestHairpinModule.UnfoldingFreeEnergy
	hairpinModule.StandbyFreeEnergyWithoutSliding = bestFreeEnergyStandby
}

func unfoldHairpin(hairpinModule HairpinModule) (HairpinModule, float64, bool) {
	stem := hairpinModule.Hairpin.Stem
	if stem.EnclosedFivePrimeIdx == -1 {
		// There are no more base pairs to unfold
		return hairpinModule, 0.0, false
	} else {
		// There is atleast one stem structure to unfold
		stemStructureToUnfold := stem.Structures[0]
		dG_unfolding := float64(-stemStructureToUnfold.Energy) / 100

		// Update stem to remove structure at base of stem
		stem.Structures = stem.Structures[1:]
		stem.ClosingFivePrimeIdx = stemStructureToUnfold.EnclosedFivePrimeIdx
		stem.ClosingThreePrimeIdx = stemStructureToUnfold.EnclosedThreePrimeIdx

		nbUnpairedFivePrimeStem := stemStructureToUnfold.EnclosedFivePrimeIdx - stemStructureToUnfold.ClosingFivePrimeIdx
		nbUnpairedThreePrimeStem := stemStructureToUnfold.ClosingThreePrimeIdx - stemStructureToUnfold.EnclosedThreePrimeIdx

		if stem.ClosingFivePrimeIdx == stem.EnclosedFivePrimeIdx {
			// stem now only consists of its closing base pair so update
			// enclosed pair to -1 to indicate no more structures in stem
			if len(stem.Structures) != 0 {
				// sanity check
				panic("hairpin stem should not contain any more structures but does")
			}
			stem.EnclosedFivePrimeIdx, stem.EnclosedThreePrimeIdx = -1, -1
		}

		// update hairpinModule's hairpin
		hairpinModule.Hairpin.Stem = stem

		// update distal and proximal sites based on nubmer of nucleotides unpaired
		hairpinModule.distalSingleStranded = addNucleotidesDistal(hairpinModule.distalSingleStranded, nbUnpairedFivePrimeStem, stemStructureToUnfold.EnclosedFivePrimeIdx-1)
		hairpinModule.proximalSingleStranded = addNucleotidesProximal(hairpinModule.proximalSingleStranded, nbUnpairedThreePrimeStem, stemStructureToUnfold.EnclosedThreePrimeIdx+1)

		return hairpinModule, dG_unfolding, true
	}
}

func addNucleotidesDistal(ssr *SingleStrandedRegion, nbNucleotides, newSSRThreePrimeIdx int) *SingleStrandedRegion {
	if ssr == nil {
		return &SingleStrandedRegion{
			ThreePrimeIdx: newSSRThreePrimeIdx,

			FivePrimeIdx: newSSRThreePrimeIdx - nbNucleotides + 1,
		}
	} else {
		ssrCopy := *ssr
		ssrCopy.ThreePrimeIdx += nbNucleotides
		return &ssrCopy
	}
}

func addNucleotidesProximal(ssr *SingleStrandedRegion, nbNucleotides, newSSRFivePrimeIdx int) *SingleStrandedRegion {
	if ssr == nil {
		return &SingleStrandedRegion{
			FivePrimeIdx: newSSRFivePrimeIdx,
			// we add -1 as the five prime idx of the ssr adds one nucleotide to the ssr
			ThreePrimeIdx: newSSRFivePrimeIdx + nbNucleotides - 1,
		}
	} else {
		ssrCopy := *ssr
		ssrCopy.FivePrimeIdx -= nbNucleotides
		return &ssrCopy
	}
}

func dG_distortion(surfaceArea int) float64 {
	if surfaceArea >= 22 {
		return 0.00
	}

	surfaceAreaFloat64 := float64(surfaceArea)
	// Original Params

	// RBS Calculator v1 - v2.0 parameters:
	c1, c2, c3 := 0.0382, -1.628, 17.3586

	// RBS Calculator v2.1 parameters:
	// c1 := 0.093  // kcal/mol/nt2
	// c2 := -1.815 // kcal/mol/nt
	// c3 := 8.857  // kcal/mol
	return c1*surfaceAreaFloat64*surfaceAreaFloat64 + c2*surfaceAreaFloat64 + c3
}

func HairpinModules(structures []interface{}) (ret []*HairpinModule) {
	for idx, structure := range structures {
		switch structure := structure.(type) {
		case *MultiLoop:
			ret = append(ret, HairpinModules(structure.Substructures)...)
		case *Hairpin:
			// find the previous structure and set it if its a single stranded region
			var prevStructurePointer *SingleStrandedRegion = nil
			if idx > 0 {
				switch prevStructure := structures[idx-1].(type) {
				case *SingleStrandedRegion:
					prevStructurePointer = prevStructure
				}
			}

			// find the next structure and set it if its a single stranded region
			var nextStructurePointer *SingleStrandedRegion = nil
			if idx < len(structures)-1 {
				switch nextStructure := structures[idx+1].(type) {
				case *SingleStrandedRegion:
					nextStructurePointer = nextStructure
				}
			}

			hairpinModule := newHairpinModule(structure, prevStructurePointer, nextStructurePointer)
			for _, prevHairPinModule := range ret {
				// add this hairpin module to the list of downstream hairpin modules of previous hairpin modules
				prevHairPinModule.downstreamHairpinModules = append(prevHairPinModule.downstreamHairpinModules, &hairpinModule)
			}
			ret = append(ret, &hairpinModule)
		}
	}

	return
}

func AddDownstreamModulesSlidingFreeEnergy(hairpinModules []*HairpinModule) []*HairpinModule {
	for _, hairpinModule := range hairpinModules {
		hairpinModule.addDgSliding()
	}
	return hairpinModules
}

func (hairpinModule *HairpinModule) addDgSliding() {
	downstreamModules := hairpinModule.downstreamHairpinModules
	var hairpinHeightDownstreamModules int
	for _, hairpinModule := range downstreamModules {
		hairpinHeightDownstreamModules += hairpinModule.HairpinHeight
	}
	slidingCoefficient := 0.2 // kcal/mol/nt
	slidingFreeEnergy := slidingCoefficient * float64(hairpinHeightDownstreamModules)
	hairpinModule.StandbyFreeEnergy = hairpinModule.StandbyFreeEnergyWithoutSliding + slidingFreeEnergy
	hairpinModule.SlidingFreeEnergy = slidingFreeEnergy
}

func (hairpinModule *HairpinModule) FivePrimeIdx() int {
	if hairpinModule.distalSingleStranded == nil {
		return hairpinModule.Hairpin.Stem.ClosingFivePrimeIdx
	} else {
		return hairpinModule.distalSingleStranded.FivePrimeIdx
	}
}

func (hairpinModule *HairpinModule) ThreePrimeIdx() int {
	if hairpinModule.proximalSingleStranded == nil {
		return hairpinModule.Hairpin.Stem.ClosingThreePrimeIdx
	} else {
		return hairpinModule.proximalSingleStranded.ThreePrimeIdx
	}
}
