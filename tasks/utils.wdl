version development

task GetGeneChrPairs {

    input {
        File featureVariantFile
    }

	command <<<

cat << EOF > script.py
import csv

with open("~{featureVariantFile}") as f, open("outputPairs.tsv", "w+") as w:
    f.readline()
    cut_lines = [(l[1], l[3]) for l in csv.reader(f)]
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
