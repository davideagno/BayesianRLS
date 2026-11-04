#!/bin/bash
# Runs all the analyses of the paper, sequentially (the fastest first), from the folder of each script.
# Usage: ./run_all.sh [folder ...]   (e.g. ./run_all.sh sec5.2_additive_models; default: all)
#
# Output goes to <folder>/output/. The log run_all.log records, for each step, the
# elapsed time and the number of system sleep events (macOS) during the step: computing
# times are reliable only for steps with 0 sleep events.

cd "$(dirname "$0")"

# On macOS, wake the machine fully and prevent idle and system sleep while running (the
# display may still turn off)
if command -v caffeinate >/dev/null && [ -z "$RUN_ALL_AWAKE" ]; then
  RUN_ALL_AWAKE=1 exec caffeinate -imsu "$0" "$@"
fi

LOG="$PWD/run_all.log"
FOLDERS="$*"

steps=(
  "secS2.1_illustration student_t_errors.R"
  "secS2.1_illustration skew_t_errors.R"
  "secS2.1_illustration posteriors_figure.R"
  "sec5.2_additive_models 1_simulate_data.R 30"
  "sec5.2_additive_models 1_simulate_data.R 100"
  "sec5.2_additive_models 2_fit_rls_student_t.R 30"
  "sec5.2_additive_models 2_fit_rls_student_t.R 100"
  "sec5.2_additive_models 3_fit_rls_gaussian.R 30"
  "sec5.2_additive_models 3_fit_rls_gaussian.R 100"
  "sec5.2_additive_models 4_fit_blm_brms.R 30"
  "sec5.2_additive_models 4_fit_blm_brms.R 100"
  "sec5.2_additive_models 5_fit_bart.R 30"
  "sec5.2_additive_models 5_fit_bart.R 100"
  "sec5.2_additive_models 6_summarize_table.R 30"
  "sec5.2_additive_models 6_summarize_table.R 100"
  "secS2.3_frequentist_am 1_fit_frequentist_am.R 30"
  "secS2.3_frequentist_am 1_fit_frequentist_am.R 100"
  "secS2.3_frequentist_am 2_summarize_table.R 30"
  "secS2.3_frequentist_am 2_summarize_table.R 100"
  "sec6_application 2_analysis.R"
  "sec5.1_uq_gaussian simulate.R 3"
  "sec5.1_uq_gaussian simulate.R 10"
  "sec5.1_uq_gaussian summarize_tables.R 3 norm"
  "sec5.1_uq_gaussian summarize_tables.R 10 norm"
  "secS2.2_uq_student_t simulate.R 3"
  "secS2.2_uq_student_t simulate.R 10"
  "sec5.1_uq_gaussian summarize_tables.R 3 t"
  "sec5.1_uq_gaussian summarize_tables.R 10 t"
)

log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

# Number of sleep events recorded by macOS between two timestamps "YYYY-MM-DD HH:MM:SS"
sleeps_between() {
  command -v pmset >/dev/null || { echo NA; return; }
  pmset -g log | awk -v a="$1" -v b="$2" \
    '/Entering Sleep/ { t = $1 " " $2; if (t >= a && t <= b) n++ } END { print n + 0 }'
}

for step in "${steps[@]}"; do
  set -- $step
  folder=$1; script=$2; shift 2
  if [ -n "$FOLDERS" ] && [[ " $FOLDERS " != *" $folder "* ]]; then continue; fi
  mkdir -p "$folder/output"
  start=$(date '+%F %T'); t0=$(date +%s)
  log "START $folder/$script $*"
  (cd "$folder" && Rscript "$script" "$@" > "output/log_${script%.R}${1:+_$1}${2:+_$2}.txt" 2>&1)
  status=$?
  end=$(date '+%F %T')
  log "END   $folder/$script $* exit=$status elapsed=$(( $(date +%s) - t0 ))s sleeps=$(sleeps_between "$start" "$end")"
done
log "ALL DONE"
