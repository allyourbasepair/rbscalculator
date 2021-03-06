package model

import (
	"fmt"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/allyourbasepair/rbscalculator/csv_helper"
	"github.com/allyourbasepair/rbscalculator/model/datasets"
)

/*****************************************************************************
Process to develop a model for the RBS Calculator

Current research provides intuition as to the process behind a ribosome binding
to a mRNA, but there is no concensus on an actual model of this interaction.

To develop a model, we need examine different primary and secondary structure
features and figure out their relationship to the end result (translation
initiation rate).

To create the model, we follow a Data -> Model Theory, Variables & Assumptions
-> Visualization & Quantitative Relationship Formalization -> Final Model
Equation workflow.


Data:
We need a high quality dataset that includes the mRNA sequence, its structure,
the organism the sequence was inserted into, the protein expressed by the CDS
and the amount of protein synthesized.
The mRNA sequence must be split into its 5' untranslated region (5' UTR) and its
protein coding sequence (CDS). This data is available in the `dataset` folder.
Please read the readme file in the folder for more information about the
datasets (./datasets/readme.md).


Model Theory, Variables & Assumptions:
After we have the data, we need to create a theoretical model of the
interactions of a ribosome with a mRNA strand. Once we have a model, we need to
write the code to come up with the variables required by the model. This part
of the process is specific to the theories people have in mind so nothing more
can be generalized about this part. If you'd like to create your own model,
create a subpackage under this subpackage (for example `model/salis_lab_v2_1`)
and describe your model and write the code to compute properties for your model
in a `properties.go` file in your created subpackage.
For an example, please have a look at the `properties.go` file in the
`salis_lab_v2_1` subpackage (./salis_lab_v2_1/properties.go).

Once the code is written to compute your desired properties, we have to run
the ground-truth dataset against your code to generate properties for all the
data points in the dataset. This package exports the function
`ComputePropertiesForDataset` which takes each data point in the
'ground-truth' dataset, computes the properties you've specified, and outputs a
csv file in a directory named `dataset_with_properties` in your subpackage.

Please see the `properties_test.go` file in the `salis_lab_v2_1` subpackage
(./salis_lab_v2_1/properties.go) to understand how to create a dataset with
your computed properties.
You can see an example of the dataset file after properties have been computed
(./salis_lab_v2_1/dataset_with_properties/train_*.csv).


Visualization & Quantitative Relationship Formalization:

After the above step, we have a dataset with properties, but we don't yet
quantitatively understand the relationship between a computed property and
the ground-truth protein levels.

Thus, to understand these relationships, we need to visualize the relationship
between a computed property and the ground-truth protein levels.

Once we can see a (linear, quadratic, cubic) relationship between a variable
and the outcome, we can figure out the quantitative relationship through
curve-fitting apps like MATLAB and add this quantitative relationship to our
model. Other quantitative analysis software like Excel can be helpful here.

Unfortunately, guidance to do this step is not yet included in this repo as
we've based our first calculator on research papers from The Salis Lab which
have included the quantitative relationship between properties and the
resulting protein levels. We hope to include this information and a workflow
for this step in the future if / when the model of the RBS calculator is
improved. PRs are very welcome.


Final Model Equation:

Although the previous step gives us the relationships between individual
variables and the outcome, it could be the case that adding all the variables
together leads to worse prediction of results.

Thus, as a final step we run regression analysis that tries out all the possible
combinations of variables in the model to see if removing some terms could
actually improve the predicted results. The least set of variables with the
highest correlation with actual results are then used in the model.

Unfortunately, as in the step above, guidance to do this is not included in
this repo as we've not had to do this yet. PRs are very welcome.


Implementation details to create a model:
To create a model, the following must be done -
* Create a subpackage in the `model` package with a `properties.go` file.
* The `properties.go` file should include a `PropertiesToCompute` variable of
  type `[]RBSPropertyFunc`.
	The value returned by each `RBSPropertyFunc` in `PropertiesToCompute` is
	added to the `Properties` map (the key in the map is the name of the
	`RBSPropertyFunc`) of the `RibosomeBindingSite` struct. Thus, the order of
	the properties in `PropertiesToCompute` is important and properties later in
	the `PropertiesToCompute` slice can use the values of properties earlier in
	the slice by accessing the relevant value using the `RBSPropertyFunc` name
	of the earlier property.
	(Since the name of the `RBSPropertyFunc` you create is likely to change, there
	is a helper func `PropertyValue` that takes a `RBSPropertyFunc`, gets its
	name as a string and then returns the value from the `Properties` field
	of a `RibosomeBindingSite` struct.)
* Create a copy of `properties_test.go` from the `salis_lab_v2_1` subpackage
	into your own subpackage. You will have to update the name of the package on
	the top of the file. Update the `propertiesToOutputToCSV` to include
	properties you created that you would like to output to the CSV.
	Finally, create the dataset with your computed properties by `cd`ing into
	the directory of your subpackage and running
	`go test -timeout 0 -run ^TestComputeProperties_forTrainDataset`.
* You will now have a csv file with your computed properties for each data
	point in the `train` dataset. At this point, you will have to curve-fit and
	quantitatively figure out the relationship between the values of each computed
	property and their relationship to the actual protein mean and protein std
	values. (As mentioned above, we don't provide guidance on this step currently.
	Excel and MATLAB are useful tools for this. PRs are very welcome.)
* Update your `properties.go` file to include the quantitative relationships
	from the previous step.
* Finally, create a function named `ComputeTranslationInitiationRate` of type
	`func(RibosomeBindingSite) float64` that returns the translation initiation
	rate for a ribosome binding site.
* Update `rbs_calculator.go` in the `rbs_calculator` package to use your model
	by updating the `rbs_model` import at the top of the file.

For an example, have a look at the `properties.go` and `properties_test.go`
files in the `salis_lab_v2_1` package.

To compute statistics for your model, read the documentation of `stats.py`.

*****************************************************************************/

// RibosomeBindingSite is a struct to represent a ribosome binding site.
type RibosomeBindingSite struct {
	// FivePrimeUTR is the untranslated region of the mRNA sequence on the five
	// prime end
	FivePrimeUTR string
	// ProteinCodingSequence is the sequence after the 5' UTR that encodes a
	// protein
	ProteinCodingSequence string
	Temperature           float64
	// RibosomalRNA is the 16S rRNA of the organism the mRNA strand is inserted
	// into
	RibosomalRNA string
	// Properties is a map that contains the computed properties
	Properties map[string]interface{}
	// OtherInformation contains information that we'd like to output to the
	// dataset with properties, but is not used to compute properties of a
	// binding site. This field is used by the `ComputePropertiesForDataset` func.
	OtherInformation []string
	// TranslationInitiationRate is the computed translation initiation rate for
	// this binding site. Please note that this field is only set during inference
	// in the `TranslationInitiationRate` func in the `rbs_calculator` package
	TranslationInitiationRate float64
}

// RBSPropertyFunc is the type of a func that can be used to compute a property
// for the `RibosomeBindingSite` struct
type RBSPropertyFunc (func(*RibosomeBindingSite) interface{})

// CleanRNA converts DNA to RNA and makes sequence upper case
func CleanRNA(rna string) string {
	rna = strings.TrimSpace(rna)
	rna = strings.ToUpper(rna)
	rna = strings.ReplaceAll(rna, "T", "U")
	return rna
}

// ComputePropertiesForDataset creates a `RibosomeBindingSite` struct for each
// data point in the input dataset, computes the desired properties (the
// `propertiesToCompute` argument), and outputs the selected properties of each
// `RibosomeBindingSite` struct (the `propertiesToOutputToCSV` argument) to a
// CSV in the `dataset_with_properties` directory in the package this function
// is called from.
//
// The input dataset may contain information that we'd like to include the
// output csv which may not be part of the `RibosomeBindingSite` struct.
// The `otherInformationMap` is an argument which specifies a map of column
// index to column header name which will be included in the outputted csv
// file.
func ComputePropertiesForDataset(datasetName string,
	propertiesToCompute, propertiesToOutputToCSV []RBSPropertyFunc,
	fivePrimeUTRColIdx, cdsColIdx, tempColIdx, ribosomalRNAColIdx int,
	columnsToKeepInOutputCSV ...datasets.DatasetColumn) {
	dataset := datasetName + ".csv"

	// Create a wait group that will only 'release' when all go subroutines
	// call `wg.Done()`
	var wg sync.WaitGroup

	// Step 1: populate a channel with `RibosomeBindingSite` structs from each row
	// of the dataset
	rbsChannel := make(chan *RibosomeBindingSite)
	// add to the wait group for the `populateMRNAChannelFromDataset` subroutine
	wg.Add(1)
	go populateRBSChannelFromDataset(dataset, rbsChannel, &wg, fivePrimeUTRColIdx, cdsColIdx, tempColIdx, ribosomalRNAColIdx, columnsToKeepInOutputCSV...)

	// Step 2: For each struct in `rbsChannel`, compute properties for it (based
	// on `propertiesToCompute`), create a `[]string` with properties that should
	// be outputted to a CSV file (based on `propertiesToOutputToCSV`) and push
	// the string slice to `csvOutputChannel`
	csvOutputChannel := make(chan []string)
	// add to the wait group for the `computeRBSProperties` subroutine
	wg.Add(1)
	go computeRBSProperties(rbsChannel, propertiesToCompute, propertiesToOutputToCSV, csvOutputChannel, &wg, columnsToKeepInOutputCSV...)

	// Step 3: Output data from `csvOutputChannel` to a CSV file
	datasetOutputFile := "./dataset_with_properties/" + datasetName + ".csv"
	// add the current unix timestamp to the file name
	datasetOutputFile = csv_helper.FileNameWithUNIXTimestamp(datasetOutputFile)
	// add to the wait group for the `csv_helper.WriteToCSV` subroutine
	wg.Add(1)
	go csv_helper.WriteToCSV(datasetOutputFile, csv_helper.CREATE, csvOutputChannel, &wg)

	// wait till all the subroutines call `wg.Done()`
	wg.Wait()
}

// populateRBSChannelFromDataset creates a `RibosomeBindingSite` struct for each
// row in the `csvFile` and adds it to `rbsChannel`
func populateRBSChannelFromDataset(csvFile string, rbsChannel chan *RibosomeBindingSite,
	wg *sync.WaitGroup, fivePrimeUTRColIdx, cdsColIdx, tempColIdx, RibosomalRNAColIdx int,
	columnsToKeepInOutputCSV ...datasets.DatasetColumn) {
	// close the rbs channel when this func returns
	defer close(rbsChannel)
	// call `wg.Done()` to decrement the wait group counter when this func returns
	defer wg.Done()

	// a channel that will hold rows of the CSV that are scanned
	csvRowsChannel := make(chan []string)
	// add to the wait group for the `csv_helper.ReadCSV` subroutine
	wg.Add(1)
	go csv_helper.ReadCSV(datasets.EmbeddedDatasetsDirectory, csvFile, csvRowsChannel, wg, 1)

	// weird Go syntax to wait for a message from a channel
	for {
		select {
		case csvRow, ok := <-csvRowsChannel:
			if !ok {
				// occurs when `csvRowsChannel` is closed by `csv_helper.ReadCSV`
				// so return this func as well as we don't have any more input
				return
			}

			if strings.TrimSpace(csvRow[fivePrimeUTRColIdx]) == "" && strings.TrimSpace(csvRow[cdsColIdx]) == "" {
				panic("dataset row doesn't contain both five prime UTR and coding sequence")
			}

			// parse temperature into `float64`
			var temp float64
			var err error
			if temp, err = strconv.ParseFloat(csvRow[tempColIdx], 64); err != nil {
				panic(err)
			}

			// get the other information which don't contribute to the
			// `RibosomeBindingSite` struct, but we want to keep in the output csv
			var otherInformationToKeep []string
			for _, datasetColumn := range columnsToKeepInOutputCSV {
				otherInformationToKeep = append(otherInformationToKeep, csvRow[datasetColumn.Idx])
			}

			// finally, create a `RibosomeBindingSite` struct with the required info
			rbs := RibosomeBindingSite{
				FivePrimeUTR:          CleanRNA(csvRow[fivePrimeUTRColIdx]),
				ProteinCodingSequence: CleanRNA(csvRow[cdsColIdx]),
				Temperature:           temp,
				RibosomalRNA:          CleanRNA(csvRow[RibosomalRNAColIdx]),
				Properties:            make(map[string]interface{}),
				OtherInformation:      otherInformationToKeep,
			}

			// add the struct to the channel for use in the `computeRBSProperties` subroutine
			rbsChannel <- &rbs

		}
	}
}

// computeRBSProperties computes properties for each `RibosomeBindingSite`
// struct in `rbsChannel`, creates a `[]string` with properties included in
// `propertiesToOutputToCSV`, and sends the `[]string` to `csvOutputChannel`
func computeRBSProperties(rbsChannel chan *RibosomeBindingSite,
	properties, propertiesToOutputToCSV []RBSPropertyFunc,
	csvOutputChannel chan []string,
	wg *sync.WaitGroup, columnsToKeepInOutputCSV ...datasets.DatasetColumn) {
	// close the csv output channel when this func returns
	defer close(csvOutputChannel)
	// call `wg.Done()` to decrement the wait group counter when this func returns
	defer wg.Done()

	var header []string
	// add the other information column headers to keep in the output csv file
	for _, datasetColumn := range columnsToKeepInOutputCSV {
		header = append(header, datasetColumn.Header)
	}
	// add the main `RibsomeBindingSite` struct fields to the header
	header = append(header, ribosomeBindingSiteStructFields()...)
	// add the names of the property functions (that need to be included in the
	// output csv) to the header
	propertiesToOutputToCSVPropertyNames := FuncNames(propertiesToOutputToCSV)
	header = append(header, propertiesToOutputToCSVPropertyNames...)

	// write the header to the csv
	csvOutputChannel <- header

	// weird Go syntax to wait for a message from a channel
	for {
		select {
		case rbs, ok := <-rbsChannel:
			if !ok {
				// occurs when `rbsChannel` is closed by `populateRBSChannelFromDataset`
				// so return this func as well as we don't have any more rbs structs
				// to process
				return
			}

			// compute properties for the rbs struct
			rbs.ComputeProperties(properties)

			// get the list of computed property values that need to be added to the
			// output csv
			propertyValues := rbs.PropertyValues(propertiesToOutputToCSVPropertyNames)

			// finally, add the rbs fields, other information and property values to
			// the output
			// note: this should be in the same order as the header row above
			var output []string = rbs.OtherInformation
			output = append(output, rbs.fieldValues()...)
			output = append(output, stringSlice(propertyValues)...)

			csvOutputChannel <- output
		}
	}
}

// rbsStructFieldsToExclude contains fields that we don't want to include
// This variable is used in the `ribosomeBindingSiteStructFields` and
// `RibosomeBindingSite.fieldValues` funcs
var rbsStructFieldsToExclude map[string]bool = map[string]bool{
	"OtherInformation":          true,
	"Properties":                true,
	"TranslationInitiationRate": true,
}

// ribosomeBindingSiteStructFields returns the fields of the
// `RibosomeBindingSite` struct excluding fields from the
// `rbsStructFieldsToExclude` map
func ribosomeBindingSiteStructFields() (fieldNames []string) {
	rbs := RibosomeBindingSite{}
	valueOfRBSStruct := reflect.ValueOf(rbs)
	typeOfRBSStruct := valueOfRBSStruct.Type()

	for i := 0; i < valueOfRBSStruct.NumField(); i++ {
		if valueOfRBSStruct.Field(i).CanInterface() {
			fieldName := typeOfRBSStruct.Field(i).Name

			// exclude fields that are in the `rbsStructFieldsToExclude` map
			if rbsStructFieldsToExclude[fieldName] {
				continue
			}
			fieldNames = append(fieldNames, typeOfRBSStruct.Field(i).Name)
		}
	}
	return
}

// fieldValues returns the values of the fields of `rbs` excluding fields from
// the `rbsStructFieldsToExclude` map
func (rbs RibosomeBindingSite) fieldValues() (fieldValues []string) {
	valueOfRBSStruct := reflect.ValueOf(rbs)
	typeOfRBSStruct := valueOfRBSStruct.Type()

	for i := 0; i < valueOfRBSStruct.NumField(); i++ {
		if valueOfRBSStruct.Field(i).CanInterface() {
			fieldName := typeOfRBSStruct.Field(i).Name

			// exclude fields that are in the `rbsStructFieldsToExclude` map
			if rbsStructFieldsToExclude[fieldName] {
				continue
			}
			fieldValue := toString(valueOfRBSStruct.Field(i).Interface())
			fieldValues = append(fieldValues, fieldValue)
		}
	}
	return
}

// FuncName returns a string of the name of the function `fn`
func FuncName(fn interface{}) string {
	// get the full name of the function with package path
	// `nameWithPackagePath` will be "<package>.<function name>"
	fnNameWithPackagePath := runtime.FuncForPC(reflect.ValueOf(fn).Pointer()).Name()

	// get only the function name
	fnNameWithPackagePathSlice := strings.Split(fnNameWithPackagePath, ".")
	functionName := fnNameWithPackagePathSlice[len(fnNameWithPackagePathSlice)-1]

	return functionName
}

// FuncNames returns the names of the function from a list of `Property`
func FuncNames(properties []RBSPropertyFunc) (names []string) {
	for _, property := range properties {
		names = append(names, FuncName(property))
	}
	return
}

// ComputeProperties computes the properties specified by the `properties`
// argument and adds the computed property value to the `Properties` field map
// of the `rbs` struct
func (rbs *RibosomeBindingSite) ComputeProperties(properties []RBSPropertyFunc) {
	for _, rbsPropertyFunc := range properties {
		propertyValue := rbsPropertyFunc(rbs)

		// add the computed property value to the `Properties` map
		propertyName := FuncName(rbsPropertyFunc)
		rbs.Properties[propertyName] = propertyValue
	}

	return
}

// PropertyValues returns the values of the properties specified by `propertyNames`
func (rbs *RibosomeBindingSite) PropertyValues(propertyNames []string) (ret []interface{}) {
	for _, propertyName := range propertyNames {
		ret = append(ret, rbs.Properties[propertyName])
	}
	return
}

// PropertyValue returns the value of the property `property`
func (rbs *RibosomeBindingSite) PropertyValue(property RBSPropertyFunc) interface{} {
	return rbs.Properties[FuncName(property)]
}

// toString returns the string value of `i`
func toString(i interface{}) string {
	return fmt.Sprint(i)
}

// stringSlice converts each variable of `interfaceSlice` to a string and
// returns the resulting string slice
func stringSlice(interfaceSlice []interface{}) (ret []string) {
	ret = make([]string, len(interfaceSlice))
	for i, v := range interfaceSlice {
		ret[i] = fmt.Sprint(v)
	}
	return
}

/*********************************************************
The following are default `RBSPropertyFunc` functions that return the values of
the fields of a `RibosomeBindingSite` struct.

These functions are used by the `rbs_calculator.TranslationInitiationRate` func
to add default properties to a `RibosomeBindingSite` struct.

`DefaultPropertiesToComputeBefore` is called before any properties are computed,
and `DefaultPropertiesToComputeAfter` is called after all properties have been
computed
*********************************************************/

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
//
// We call this before any properties have been computed as the five prime
// UTR could be updated by a property which would give us an incorrect start
// position
func StartPosition(rbs *RibosomeBindingSite) interface{} {
	return len(rbs.FivePrimeUTR)
}

// DefaultPropertiesToComputeBefore specifies the properties to compute
// before any of the `PropertiesToCompute` properties of the RBS calc model are
// computed
var DefaultPropertiesToComputeBefore []RBSPropertyFunc = []RBSPropertyFunc{
	ensureSaneInput,
	StartPosition,
}

// FivePrimeUTR returns the `FivePrimeUTR` field of the `rbs` struct
func FivePrimeUTR(rbs *RibosomeBindingSite) interface{} {
	return rbs.FivePrimeUTR
}

// ProteinCodingSequence returns the `ProteinCodingSequence` field of the `rbs`
// struct
func ProteinCodingSequence(rbs *RibosomeBindingSite) interface{} {
	return rbs.ProteinCodingSequence
}

// TranslationInitiationRate returns the `TranslationInitiationRate` field of
// the `rbs` struct
func TranslationInitiationRate(rbs *RibosomeBindingSite) interface{} {
	return rbs.TranslationInitiationRate
}

// DefaultPropertiesToComputeAfter specifies the properties to compute
// after all of the `PropertiesToCompute` properties of the RBS calc model have
// been computed
var DefaultPropertiesToComputeAfter []RBSPropertyFunc = []RBSPropertyFunc{
	FivePrimeUTR,
	ProteinCodingSequence,
	TranslationInitiationRate,
}

/******************************************************************************
End of default properties
******************************************************************************/
