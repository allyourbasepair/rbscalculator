// +build gen_csv

// Since we don't want the tests below to run when `go test ./...` is used, the
// above build contrain flag ensures the tests below are only run when we pass
// the `gen_csv` tag to go test with the tags argument
// (`go test ./... -tags gen_csv`).

package datasets

func ExampleAdd16SrRNA() {
	datasetName := "train"
	organismColIdx := 2
	add16SrRNA(datasetName, organismColIdx)

	// Output:
	// confirm results and rename "./with_16s_rRNA/<datasetName>_<timestamp>.csv" to "<datasetName>.csv"
}
