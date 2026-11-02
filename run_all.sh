#!/bin/bash
# Runs all the analyses of the paper, sequentially, from the folder of each script.
# Usage: ./run_all.sh [folder ...]   (e.g. ./run_all.sh sec5.2_additive_models; default: all)
#
# Output goes to <folder>/output/. The log run_all.log records, for each step, the
# elapsed time and the number of system sleep events (macOS) during the step: computing
# times are reliable only for steps with 0 sleep events.

cd "$(dirname "$0")"

# On macOS, keep the machine fully awake (display, idle and system sleep) while running
if command -v caffeinate >/dev/null && [ -z "$RUN_ALL_AWAKE" ]; then
  RUN_ALL_AWAKE=1 exec caffeinate -dimsu "$0" "$@"
fi

LOG="$PWD/run_all.log"
FOLDERS="$*"

steps=(
  "secS2.1_illustration student_t.R"
  "secS2.1_illustration skew_t.R"
  "sec5.1_uq_gaussian sim_normal.R 3"
  "sec5.1_uq_gaussian sim_normal.R 10"
  "sec5.1_uq_gaussian results.R 3 norm"
  "sec5.1_uq_gaussian results.R 10 norm"
  "secS2.2_uq_student_t sim_t.R 3"
  "secS2.2_uq_student_t sim_t.R 10"
  "sec5.1_uq_gaussian results.R 3 t"
  "sec5.1_uq_gaussian results.R 10 t"
  "sec5.2_additive_models 1_sim_data.R 30"
  "sec5.2_additive_models 1_sim_data.R 100"
  "sec5.2_additive_models 2_rls.R 30"
  "sec5.2_additive_models 2_rls.R 100"
  "sec5.2_additive_models 2_rls_norm.R 30"
  "sec5.2_additive_models 2_rls_norm.R 100"
  "sec5.2_additive_models 3_trad.R 30"
  "sec5.2_additive_models 3_trad.R 100"
  "sec5.2_additive_models 4_bart.R 30"
  "sec5.2_additive_models 4_bart.R 100"
  "sec5.2_additive_models 5_results.R 30"
  "sec5.2_additive_models 5_results.R 100"
  "secS2.3_frequentist_am freq_AM.R 30"
  "secS2.3_frequentist_am freq_AM.R 100"
  "secS2.3_frequentist_am freq_AM_results.R 30"
  "secS2.3_frequentist_am freq_AM_results.R 100"
  "sec6_application analysis.R"
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
