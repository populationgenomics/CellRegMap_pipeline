version development

# additional useful tasks (not necessarily linked to CellRegMap)
# this task take a csv file and extracts pairs of values based on specific column names

task CsvPairExtractor {

    input {
        File csvFile
        Array[String]+ columnsToSelect # chrom, gene
    }

	command <<<

cat << EOF > script.py
import csv

columns_to_select = '~{sep(",", columnsToSelect)}'.split(",")
# we expect output PAIRS, so make sure there are 2 columns to select
assert len(columns_to_select) == 2
with open("~{csvFile}") as f, open("outputPairs.tsv", "w+") as w:
    reader = csv.reader(f)

    # get headers list to get a header to indices
    headers = next(reader)
    # this will fail if column isn't in the header
    indices_to_select = [headers.index(h) for h in columns_to_select]

    cut_lines = [[l[i] for i in indices_to_select] for l in reader]
    unique_lines = set('\t'.join(l) for l in cut_lines)
    if unique_lines:
    # this should leave an empty fails
        w.write("\n".join(sorted(unique_lines)))

EOF

python script.py

	>>>

    runtime {
        memory: "2GB"
        docker: "python:3.10"
    }

    output {
        Array[Array[String]] outputPairs = read_tsv("outputPairs.tsv")
    }
}
