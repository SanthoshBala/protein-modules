package org.hupo.psi.mi.psicquic.indexing.batch;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.core.Job;
import org.springframework.batch.core.JobExecution;
import org.springframework.batch.core.JobParametersBuilder;
import org.springframework.batch.core.JobParametersInvalidException;
import org.springframework.batch.core.launch.*;
import org.springframework.batch.core.repository.JobExecutionAlreadyRunningException;
import org.springframework.batch.core.repository.JobInstanceAlreadyCompleteException;
import org.springframework.batch.core.repository.JobRepository;
import org.springframework.batch.core.repository.JobRestartException;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationContext;
import org.springframework.context.support.ClassPathXmlApplicationContext;

import javax.annotation.Resource;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * Main class to index mitab in solr
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/05/12</pre>
 */

public class SolrMitabIndexer {

    private static final Log log = LogFactory.getLog(SolrMitabIndexer.class);

    @Resource(name = "batchJobLauncher")
    private JobLauncher jobLauncher;

    @Resource(name = "jobRepository")
    private JobRepository jobRepository;

    @Resource(name = "jobOperator")
    private JobOperator jobOperator;

    @Autowired
    private ApplicationContext applicationContext;

    private String indexingId;

    public SolrMitabIndexer() {
    }

    public static void main(String[] args) throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, NoSuchJobExecutionException, NoSuchJobException, NoSuchJobInstanceException {

        // loads the spring context defining beans and jobs
        ApplicationContext context = new ClassPathXmlApplicationContext(
                new String[] {"/META-INF/psicquic-spring.xml", "/META-INF/jobs/psicquic-indexing-spring.xml"});

        SolrMitabIndexer rm = (SolrMitabIndexer)
                context.getBean("solrMitabIndexer");
        
        if ( args.length == 1){
            rm.setIndexingId(args[0]);
            rm.resumeIndexing();
        }
        else {
            rm.startIndexing();
        }
    }

    public void resumeIndexing() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, NoSuchJobExecutionException, JobRestartException, NoSuchJobException, NoSuchJobInstanceException {
        resumeJob("mitabIndexJob");
    }

    public void resumeJob(String jobName) throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, NoSuchJobExecutionException, JobRestartException, NoSuchJobException, NoSuchJobInstanceException {
        Long executionId = findJobId(indexingId, jobName);

        if (executionId == null) {
            throw new IllegalStateException("Indexing Id not found: "+indexingId);
        }

        jobOperator.restart(executionId);
    }

    private Long findJobId(String indexingId, String jobName) throws NoSuchJobException, NoSuchJobExecutionException, NoSuchJobInstanceException {
        final List<Long> jobIds = jobOperator.getJobInstances(jobName, 0, Integer.MAX_VALUE);

        for (Long jobId : jobIds) {
            for (Long executionId : jobOperator.getExecutions(jobId)) {
                final String params = jobOperator.getParameters(executionId);

                if (params.contains(indexingId)) {
                    return executionId;
                }
            }
        }

        return null;
    }

    public void startIndexing() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException {
        startJob("mitabIndexJob");
    }

    public void startJob(String jobName) throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException {
        String indexingId = "psicquic_index_"+System.currentTimeMillis();
        setIndexingId(indexingId);

        if (log.isInfoEnabled()) log.info("Starting indexing : "+indexingId);

        FileWriter indexingWriter = null;
        try {
            String fileName = "target/"+jobName+"_"+System.currentTimeMillis();
            if (log.isInfoEnabled()) log.info("The indexing id is stored in : "+fileName +". This indexing id is necessary for restarting an indexing job from where it failed.");

            indexingWriter = new FileWriter(fileName);

            indexingWriter.write(indexingId);
        } catch (IOException e) {
            log.error("Impossible to create the file containing the indexing id for this job.", e);
        }
        finally {
            try {
                indexingWriter.close();
            } catch (IOException e) {
                log.error("Impossible to close the file writer writing the indexing id for this job.", e);
            }
        }

        runJob(jobName, indexingId);
    }

    protected JobExecution runJob(String jobName, String indexingId) throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException {
        if (log.isInfoEnabled()) log.info("Starting job: "+jobName);
        Job job = (Job) applicationContext.getBean(jobName);

        JobParametersBuilder jobParamBuilder = new JobParametersBuilder();
        jobParamBuilder.addString("indexingId", indexingId).toJobParameters();

        log.info("starting job " + indexingId);

        return jobLauncher.run(job, jobParamBuilder.toJobParameters());
    }

    public void setIndexingId(String indexingId) {
        this.indexingId = indexingId;
    }

    public String getIndexingId() {
        return indexingId;
    }
}
