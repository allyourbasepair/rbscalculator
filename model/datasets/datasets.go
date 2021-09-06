package datasets

import (
	"embed"
	"sync"

	"github.com/allyourbasepair/rbscalculator/csv_helper"
)

/***********************************************
Functions for modifying datasets in the `datasets` directory.

Currently only includes a func `add16SrRNA` to add the 16S
rRNA of an organism to a dataset.
***********************************************/

//go:embed *
var EmbeddedDatasetsDirectory embed.FS

// DatasetColumn represents a column in a dataset. This is used to keep dataset
// columns in the output CSV when generating a CSV with computed properties.
type DatasetColumn struct {
	Header string
	Idx    int
}

// add16SrRNA adds the 16s rRNA (based on the
// `rbs_calculator.Organism16SrRNAMap` map) to a dataset
func add16SrRNA(datasetName string, organismColIdx int) {
	dataset := datasetName + ".csv"

	// Create a wait group that will only 'release' when all go subroutines
	// call `wg.Done()`
	var wg sync.WaitGroup

	nbHeaderRows := 1
	headerRow := csv_helper.ReadHeader(EmbeddedDatasetsDirectory, dataset, nbHeaderRows)[0]
	headerRow = append(headerRow, "16S rRNA")

	// Step 1: populate a channel with each row of the dataset
	csvRowsChannel := make(chan []string)
	// add to the wait group for the `csv_helper.ReadCSV` subroutine
	wg.Add(1)
	go csv_helper.ReadCSV(EmbeddedDatasetsDirectory, dataset, csvRowsChannel, &wg, nbHeaderRows)

	// Step 2: For each row in the dataset, add the 16s rRNA of the organism
	csvOutputChannel := make(chan []string, 1)
	// write the header row to file
	csvOutputChannel <- headerRow
	// add to the wait group for the `doAdd16SrRNA` subroutine
	wg.Add(1)
	go doAdd16SrRNA(csvRowsChannel, organismColIdx, csvOutputChannel, &wg)

	// Step 3: Output data from `csvOutputChannel` to a CSV file
	datasetOutputFile := "./with_16s_rRNA/" + datasetName + ".csv"
	// add the current unix timestamp to the file name
	datasetOutputFile = csv_helper.FileNameWithUNIXTimestamp(datasetOutputFile)
	// add to the wait group for the `csv_helper.WriteToCSV` subroutine
	wg.Add(1)
	go csv_helper.WriteToCSV(datasetOutputFile, csv_helper.CREATE, csvOutputChannel, &wg)

	// wait till all the subroutines call `wg.Done()`
	wg.Wait()
}

func doAdd16SrRNA(csvInputChannel chan []string, organismColIdx int, csvOutputChannel chan []string, wg *sync.WaitGroup) {
	// close the csv output channel when this func returns
	defer close(csvOutputChannel)
	// call `wg.Done()` to decrement the wait group counter when this func returns
	defer wg.Done()

	// weird Go syntax to wait for a message from a channel
	for {
		select {
		case csvRow, ok := <-csvInputChannel:
			if !ok {
				return
			}

			// get the organism for the current dataset row
			// organism := csvRow[organismColIdx]

			// if val, ok := rbs_calculator.Organism16SrRNAMap[organism]; ok {
			// 	csvRow = append(csvRow, val)
			// } else {
			// 	panic(fmt.Errorf("do not have 16S rRNA for organism: %v", organism))
			// }

			csvOutputChannel <- csvRow
		}
	}
}
