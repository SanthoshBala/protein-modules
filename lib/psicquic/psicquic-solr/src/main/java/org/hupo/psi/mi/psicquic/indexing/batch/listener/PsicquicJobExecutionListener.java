package org.hupo.psi.mi.psicquic.indexing.batch.listener;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.core.JobExecution;
import org.springframework.batch.core.JobExecutionListener;

/**
 * Listener of a PSICQUIC Job
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/05/12</pre>
 */

public class PsicquicJobExecutionListener implements JobExecutionListener {

    private static final Log log = LogFactory.getLog(SimpleChunkListener.class);

    public void beforeJob(JobExecution jobExecution) {
        log.info("\nJOB STARTED: " + jobExecution + "\n");
    }

    public void afterJob(JobExecution jobExecution) {
        log.info("\nJOB FINISHED: "+jobExecution+"\n");
    }
}
