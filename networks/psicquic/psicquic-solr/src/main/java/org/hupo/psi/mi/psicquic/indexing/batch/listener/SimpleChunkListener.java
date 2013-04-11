package org.hupo.psi.mi.psicquic.indexing.batch.listener;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.core.ChunkListener;
import org.springframework.batch.core.ExitStatus;
import org.springframework.batch.core.StepExecution;
import org.springframework.batch.core.StepExecutionListener;

/**
 * Listener of a job step that will print the number of lines that have been read
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/05/12</pre>
 */

public class SimpleChunkListener implements StepExecutionListener, ChunkListener {

    private static final Log log = LogFactory.getLog(SimpleChunkListener.class);

    private StepExecution stepExecution;
    private long startTime;

    public ExitStatus afterStep(StepExecution stepExecution) {
        final int readCount = stepExecution.getReadCount();

        return null;
    }

    public void beforeStep(StepExecution stepExecution) {
        this.stepExecution = stepExecution;

        this.startTime = System.currentTimeMillis();
    }

    public void beforeChunk() {
    }

    public void afterChunk() {
        if (!log.isInfoEnabled()) {
            return;
        }

        final int readCount = stepExecution.getReadCount();

        log.info("Number of lines read : " + readCount);
    }
}
