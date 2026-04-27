version 1.1

workflow main {
	input {
		File msFile
		String msSheet
		String msSite
		String outDir
		String date
		String bait
		String cellType
		String samples
		String controls
		String baitInWeb
	}
	
	call run_genoppi_saintexp {
		input:
			msFile = msFile,
			msSheet = msSheet,
			msSite = msSite,
			outDir = outDir,
			date = date,
			bait = bait,
			cellType = cellType,
			samples = samples,
			controls = controls,
			baitInWeb = baitInWeb
	}
	
	output {
		File outFile = run_genoppi_saintexp.outTarball
	}
}

task run_genoppi_saintexp {
	input {
		File msFile
		String msSheet
		String msSite
		String outDir
		String date
		String bait
		String cellType
		String samples
		String controls
		String baitInWeb
	}

	command <<<
	
		Rscript /usr/local/src/runGenoppi.r \
		"~{msFile}" \
		"~{msSheet}" \
		"~{msSite}" \
		"~{outDir}" \
		"~{date}" \
		"~{bait}" \
		"~{cellType}" \
		"~{samples}" \
		"~{controls}" \
		"~{baitInWeb}"
		
		Rscript /usr/local/src/runSAINTexpressInt.r \
		"~{msFile}" \
		"~{msSheet}" \
		"~{msSite}" \
		"~{outDir}" \
		"~{date}" \
		"~{bait}" \
		"~{cellType}" \
		"~{samples}" \
		"~{controls}" \
		"~{baitInWeb}"
		
		Rscript /usr/local/src/runSAINTexpressSpc.r \
		"~{msFile}" \
		"~{msSheet}" \
		"~{msSite}" \
		"~{outDir}" \
		"~{date}" \
		"~{bait}" \
		"~{cellType}" \
		"~{samples}" \
		"~{controls}" \
		"~{baitInWeb}"
		
		tar -czvf "~{outDir}.tar.gz" "~{outDir}"
	>>>

	output {		
		File outTarball = "~{outDir}.tar.gz"
	}
	
	runtime {
		docker: "us-central1-docker.pkg.dev/lage-genoppi/genoppi/genoppi-saintexp"
		memory: "2 GB"
		cpu: 1
		preemptible: 1
	}
}

