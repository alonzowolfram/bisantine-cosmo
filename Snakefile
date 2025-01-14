## Snakemake - bisantine-cosmo
##
## @alonzowolfram
##

# --- Necessary Python packages --- #
from datetime import datetime
import sys 
import os
import filecmp
import shutil

# --- Importing configuration files --- #
# https://stackoverflow.com/questions/67108673/accessing-the-path-of-the-configfile-within-snakefile
args = sys.argv
CONFIG_PATH = args[args.index("--configfiles") + 1]
configfile: CONFIG_PATH

# --- Setting variables --- #
# Output path
def generateOutputPath(previous_run_out_dir, output_path, project_name, run_name, now):
    # Set up the project_name and run_name strings.
    if project_name is None or project_name=="":
        project_name_string = ""
    else:
        project_name_string = project_name + "_"
    if run_name is None or run_name=="":
        run_name_string = ""
    else:
        run_name_string = run_name + "_"

    if previous_run_out_dir is None or previous_run_out_dir == "":
        if output_path is None or output_path == "":
            output_path = "out/" + project_name_string + run_name_string + str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "-" + str(now.hour) + "-" + str(now.minute) + "-" + str(now.second) + "/"
        else:
            output_path = output_path
    else:
        output_path = previous_run_out_dir + "/"

    return output_path

now = datetime.now()
OUTPUT_PATH = generateOutputPath(config["data"]["previous_run_out_dir"], config["output"]["output_dir"], config["project"]["meta"]["project_name"], config["project"]["meta"]["run_name"], now)

# Nextflow vs Snakemake stuff
WORKFLOW_SYSTEM = "Snakemake"
PROJECT_DIRECTORY = ""

# --- Rules --- # 
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#onstart-onsuccess-and-onerror-handlers
onsuccess:
    workflow_system = WORKFLOW_SYSTEM,
    config_path = CONFIG_PATH,
    output_path = OUTPUT_PATH,
        
    # Ensure the output directory exists.
    os.makedirs(output_path[0] + "config/", exist_ok=True)

    # Get the current timestamp.
    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    # Export the current conda environment.
    os.system("conda list --export > " + output_path[0] + "config/conda_" + timestamp + ".env")

    # Extract the base name of the input config file (without extension.)
    config_name, config_ext = os.path.splitext(os.path.basename(config_path[0]))

    # Path to the new file.
    new_file = os.path.join(output_path[0] + "config/", f"{config_name}_{timestamp}{config_ext}")

    # Find the latest file in the `config` directory that matches this config base name
    existing_files = sorted(
        [f for f in os.listdir(output_path[0] + "config/") if f.startswith(config_name) and f.endswith(config_ext)]
    )
    latest_file = os.path.join(output_path[0] + "config/", existing_files[-1]) if existing_files else None

    # Compare the current YAML file with the latest file in the folder.
    if not latest_file or not filecmp.cmp(config_path[0], latest_file, shallow=False):
        # Copy the file if there are differences or no previous file exists.
        shutil.copy(config_path[0], new_file)
        print(f"Copied {config_path[0]} to {new_file}")
    else:
        print("No changes detected in configuration YAML file. Skipping copy.")

# rule export_env:
#     input:
#         CONFIG_PATH # lambda wildcards: configfile  # Dynamically use the configfile provided by the user.
#     output:
#          temp(OUTPUT_PATH + "config/.export_env_done")
#         # directory(OUTPUT_PATH  + "config"),
#         # env_file = OUTPUT_PATH + "config/conda.env"
#     params:
#         workflow_system = WORKFLOW_SYSTEM,
#         config_path = CONFIG_PATH,
#         output_path = OUTPUT_PATH,
#         config_file = CONFIG_PATH + "config.yaml"
#     priority: 100  # Ensure this runs before other rules.
#     run:
#         # Ensure the output directory exists.
#         os.makedirs(output[0], exist_ok=True)

#         # Export the current conda environment.
#         os.system("conda list --export > " + output_path + "config/conda.env")

#         # Extract the base name of the input config file (without extension.)
#         config_name, config_ext = os.path.splitext(os.path.basename(input[0]))

#         # Get the current timestamp.
#         timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

#         # Path to the new file.
#         new_file = os.path.join(output[0], f"{config_name}_{timestamp}{config_ext}")

#         # Find the latest file in the `config` directory that matches this config base name
#         existing_files = sorted(
#             [f for f in os.listdir(output[0]) if f.startswith(config_name) and f.endswith(config_ext)]
#         )
#         latest_file = os.path.join(output[0], existing_files[-1]) if existing_files else None

#         # Compare the current YAML file with the latest file in the folder.
#         if not latest_file or not filecmp.cmp(input[0], latest_file, shallow=False):
#             # Copy the file if there are differences or no previous file exists.
#             shutil.copy(input[0], new_file)
#             print(f"Copied {input[0]} to {new_file}")
#         else:
#             print("No changes detected in config.yaml. Skipping copy.")

#         # Mark the rule as done for this session.
#         with open(output[0], "w") as f:
#             f.write("Done.")

# rule make_report:
#     input: 
#         R_file_qc_metrics = OUTPUT_PATH + "Rdata/qc_metrics.rds"
#     output:
#         # test_file = OUTPUT_PATH + "pubs/test.rds"
#         report_file = OUTPUT_PATH + "pubs/report.html"
#     params:
#         workflow_system = WORKFLOW_SYSTEM,
#         script = "src/make_report.R",
#         output_path = OUTPUT_PATH,
#         current_module = "make_report",
#         config_path = CONFIG_PATH,
#         project_directory = PROJECT_DIRECTORY,
#         report_template_path = "ext/report_template.Rmd",
#         report_output_path = OUTPUT_PATH + "pubs/report.html"
#     log:
#         out = OUTPUT_PATH + "logs/make_report.out",
#         err = OUTPUT_PATH + "logs/make_report.err" 
#     shell:
#         """
#         mkdir -p {params.output_path}/pubs
#         Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file_qc_metrics} {params.project_directory} 1> {log.out} 2> {log.err}
#         """

rule normalization:
    input:
        R_file = OUTPUT_PATH + "Rdata/seuratObject_filtered.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/seuratObject_normalized.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/normalization.R",
        output_path = OUTPUT_PATH,
        current_module = "normalization",
        config_path = CONFIG_PATH,
        project_directory = PROJECT_DIRECTORY
    log:
        out = OUTPUT_PATH + "logs/normalization.out",
        err = OUTPUT_PATH + "logs/normalization.err" 
    shell:
        """
        Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} {params.project_directory} 1> {log.out} 2> {log.err}
        """

rule qc:
    input:
        R_file = OUTPUT_PATH + "Rdata/seuratObject_raw.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/seuratObject_filtered.rds",
        R_file_neg_probes = OUTPUT_PATH + "Rdata/neg-probes_filtered.rds",
        R_file_qc_metrics = OUTPUT_PATH + "Rdata/qc_metrics.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/qc.R",
        output_path = OUTPUT_PATH,
        current_module = "qc",
        config_path = CONFIG_PATH,
        project_directory = PROJECT_DIRECTORY
    log:
        out = OUTPUT_PATH + "logs/qc.out",
        err = OUTPUT_PATH + "logs/qc.err" 
    shell:
        """
        Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} {params.project_directory} 1> {log.out} 2> {log.err}
        """

rule data_import_cleaning:
    output:
        R_file = OUTPUT_PATH + "Rdata/seuratObject_raw.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/data_import_cleaning.R",
        output_path = OUTPUT_PATH,
        current_module = "data_import_cleaning",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/data-import-cleaning.out",
        err = OUTPUT_PATH + "logs/data-import-cleaning.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} 1> {log.out} 2> {log.err}"

rule test:
    output:
        R_file = OUTPUT_PATH + "Rdata/test.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/test.R",
        output_path = OUTPUT_PATH,
        current_module = "test",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/test.out",
        err = OUTPUT_PATH + "logs/test.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} 1> {log.out} 2> {log.err}"