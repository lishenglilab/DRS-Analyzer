def shellQuote(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process ADD_DECOY {
    tag "add_decoy"
    publishDir "${params.outdir}/protein_ref", mode: 'copy'
    stageInMode 'copy'

    input:
    path 'target.fasta', arity: '1'
    val philosopher_path

    output:
    path "*-decoys-*.fas", arity: '1', emit: decoy_db

    script:
    def philosopher = shellQuote(philosopher_path)
    """
    mkdir decoy_workspace
    cd decoy_workspace
    test -x ${philosopher} || { echo 'Philosopher is not executable in the task environment' >&2; exit 1; }
    cp -- ../target.fasta .
    ${philosopher} workspace --init --nocheck --temp ..
    ${philosopher} database --custom ./target.fasta
    ${philosopher} workspace --clean --nocheck

    shopt -s nullglob
    candidates=(*-decoys-*.fas)
    if (( \${#candidates[@]} != 1 )); then
        echo 'Philosopher must produce exactly one target-decoy FASTA' >&2
        exit 1
    fi
    if ! awk '/^>rev_.+/ { decoy=1; next } /^>[^[:space:]]/ { target=1 } END { exit !(decoy && target) }' "./\${candidates[0]}"; then
        echo 'Philosopher FASTA must contain target and rev_ decoy headers' >&2
        exit 1
    fi
    cp -- "\${candidates[0]}" ..
    """
}

process FRAGPIPE {
    tag "fragpipe"
    publishDir "${params.outdir}/protein", mode: 'copy'
    stageInMode 'copy'

    cpus 10
    memory '32 GB'
    time '48 h'

    input:
    path 'decoy.fasta'
    path 'input.workflow'
    path spectra, stageAs: 'spectra/run????/*', arity: '1..*'
    val manifest_fields
    val runtime

    output:
    path 'fragpipe_output', emit: fragpipe_output

    script:
    def manifestRows = spectra.withIndex().collect { spectrum, index ->
        def inputPath = shellQuote(spectrum)
        def fields = shellQuote(manifest_fields[index].join('\t'))
        """
        test -r "\$task_dir"/${inputPath} || { echo 'Staged mzML is unreadable in the task environment' >&2; exit 1; }
        printf '%s/%s\\t%s\\n' "\$task_dir" ${inputPath} ${fields} >> manifest.fp-manifest
        """
    }.join('\n')
    def fragpipe = shellQuote(runtime.launcher)
    def flags = runtime.args.collect { shellQuote(it) }.join(' ')
    def ramFlag = ''
    if (runtime.cli == '23.1') {
        def taskGb = task.memory?.toGiga()
        if (!taskGb || taskGb < 1 || taskGb > Integer.MAX_VALUE) {
            error 'ERROR: FRAGPIPE task memory must be at least 1 GiB and fit in a positive integer'
        }
        if (runtime.ramGb != null && runtime.ramGb > taskGb) {
            error "ERROR: --fragpipe_ram_gb ${runtime.ramGb} exceeds the FRAGPIPE task memory (${taskGb} GiB)"
        }
        ramFlag = "--ram ${runtime.ramGb ?: taskGb}"
    }
    def dryRunFlag = runtime.dryRun ? '--dry-run' : ''
    def runtimeChecks = runtime.executable.collect { path ->
        "test -x ${shellQuote(path)} || { printf 'Tool is not executable in the task environment: %s\\n' ${shellQuote(path)} >&2; exit 1; }"
    }.join('\n') + '\n' + runtime.readable.collect { path ->
        "test -r ${shellQuote(path)} || { printf 'Tool file is unreadable in the task environment: %s\\n' ${shellQuote(path)} >&2; exit 1; }"
    }.join('\n')
    def outputChecks = runtime.dryRun ? '''
    if grep -Eiq '(^|[[:space:]])(ERROR|FATAL)([[:space:]:]|$)|Permission denied|AccessDeniedException|not executable' fragpipe.log; then
        cat fragpipe.log >&2
        echo 'FragPipe dry-run reported an error' >&2
        exit 1
    fi
    if ! grep -Fqx 'FragPipe version 23.1' fragpipe.log ||
       ! grep -Fq "It's a dry-run, not running the commands." fragpipe.log; then
        echo 'FragPipe did not confirm a 23.1 dry-run' >&2
        exit 1
    fi
    for report in fragpipe.job fragpipe.workflow fragpipe-files.fp-manifest; do
        if [[ ! -s "$report" ]]; then
            printf 'FragPipe dry-run did not create a nonempty %s\\n' "$report" >&2
            exit 1
        fi
    done
    ''' : '''
    if [[ -z "$(find . -type f \\( -name '*.tsv' -o -name '*.csv' \\) -size +0c -print -quit)" ]]; then
        echo 'FragPipe produced no nonempty TSV or CSV output' >&2
        exit 1
    fi
    '''
    """
    ${runtimeChecks}
    task_dir=\$(pwd -P)
    mkdir fragpipe_workspace
    cd fragpipe_workspace
    : > manifest.fp-manifest
    ${manifestRows}

    export DRS_DECOY_DB_PATH="\$task_dir/decoy.fasta"
    if ! awk '
        /^database[.]db-path=/ { db++; print "database.db-path=" ENVIRON["DRS_DECOY_DB_PATH"]; next }
        /^database[.]decoy-tag=/ { tags++; if (\$0 ~ /^database[.]decoy-tag=rev_\\r?\$/) rev++ }
        { print }
        END { exit !(db == 1 && tags == 1 && rev == 1) }
    ' ../input.workflow > workflow.workflow; then
        echo 'FragPipe workflow needs exactly one database.db-path and database.decoy-tag=rev_' >&2
        exit 1
    fi

    ${fragpipe} --headless \\
        --workflow workflow.workflow \\
        --manifest manifest.fp-manifest \\
        --workdir . \\
        ${flags} ${dryRunFlag} ${ramFlag} \\
        --threads ${task.cpus} > fragpipe.log 2>&1 || {
            status=\$?
            cat fragpipe.log >&2
            exit "\$status"
        }

    test -s fragpipe.log || { echo 'FragPipe did not write a log' >&2; exit 1; }
    ${outputChecks}
    cd ..
    mv -- fragpipe_workspace fragpipe_output
    """
}
