package model

import (
	"embed"
)

func ExampleAdd16SrRNA() {
	var datasetsDirectory = "datasets"
	var embeddedDatasetsDirectory embed.FS
	datasetName := "train"
	dataset := datasetsDirectory + "/" + datasetName + ".csv"
	outputCSVFile := datasetsDirectory + "/" + datasetName + "_with_16S.csv"
	organismColIdx, rrnaColIdx := 2, 9
	add16SrRNA(embeddedDatasetsDirectory, dataset, outputCSVFile, organismColIdx, rrnaColIdx)

	// Output:
	// confirm results and rename "<datasetName>_with_16S.csv" to "<datasetName>.csv"
}
