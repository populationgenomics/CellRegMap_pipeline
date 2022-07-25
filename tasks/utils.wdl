version development

task GetGeneChrPairs {

    input {
        File featureVariantFile
        Array[String]+ columnsToSelect
    }

	command <<<

cat << EOF > script.py
import csv

columns_to_select = '~{sep(",", columnsToSelect)}'.split(",")
# we expect output PAIRS, so make sure there are 2 columns to select
assert len(columns_to_select) == 2
with open("~{featureVariantFile}") as f, open("outputPairs.tsv", "w+") as w:
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
        memory: "2G"
		container: "python:3.10"
    }

    output {
        # [["GeneName1", "chr1"], ["GeneName2","chr2"]]
        Array[Array[String]] output_pairs = read_tsv("outputPairs.tsv")
    }
}
