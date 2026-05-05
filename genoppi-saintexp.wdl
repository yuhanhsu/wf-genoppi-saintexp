version 1.0

workflow main {
	input {
		# docker image with Genoppi, SAINTexpress, and gcloud CLI
		String dockerImage = "us-central1-docker.pkg.dev/lage-genoppi/genoppi/genoppi-saintexp"
		String dockerTag = "2026.05.05"
		
		# bucket to store output files
		String destination
		
		# input arguments for R scripts
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
			dockerImage = dockerImage,
			dockerTag = dockerTag,
			destination = destination,
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
		# bucket with output files
		String outLink = run_genoppi_saintexp.outLink
	}
}

task run_genoppi_saintexp {
	input {
		String dockerImage
		String dockerTag
		String destination
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
		echo "### run Genoppi R script"	
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
		
		echo "### run SAINTexpress-int R script"
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
		
		echo "### run SAINTexpress-spc R script"
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
		
		echo "### upload output directory to destination bucket"
		gcloud storage cp -r "~{outDir}" "~{destination}/~{ourDir}"

		echo "~{destination}/~{ourDir}" > outLink.txt
	>>>

	output {
		String outLink = read_string("outLink.txt")
	}
	
	runtime {
		docker: "~{dockerImage}:~{dockerTag}"
		memory: "2 GB"
		cpu: 1
		preemptible: 1
	}
}

