from sch_simulation.amis_integration import amis_integration
from sch_simulation.helsim_FUNC_KK.file_parsing import (
    parse_coverage_input,
)

fp = amis_integration.FixedParameters(
    number_hosts=10,
    coverage_file_name="mansoni_coverage_scenario_0.xlsx",
    demography_name="UgandaRural",
    survey_type="KK2",
    parameter_file_name="mansoni_params.txt",
    coverage_text_file_storage_name="Man_MDA_vacc.txt",
    min_multiplier=5,
)

parse_coverage_input(
    fp.coverage_file_name,
    fp.coverage_text_file_storage_name,
)

amis_integration.run_and_extract_results(
    parameter_set=(3, 0.04),
    seed=1,
    fixed_parameters=fp,
    year_indices=[0, 23],
)
