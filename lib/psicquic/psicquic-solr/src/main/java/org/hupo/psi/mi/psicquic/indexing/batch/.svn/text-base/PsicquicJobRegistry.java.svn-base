package org.hupo.psi.mi.psicquic.indexing.batch;

import org.springframework.batch.core.Job;
import org.springframework.batch.core.configuration.DuplicateJobException;
import org.springframework.batch.core.configuration.JobFactory;
import org.springframework.batch.core.configuration.ListableJobLocator;
import org.springframework.batch.core.launch.NoSuchJobException;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationContext;

import java.util.Arrays;
import java.util.Collection;

/**
 * The spring job registry for PSICQUIC
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>29/05/12</pre>
 */

public class PsicquicJobRegistry implements ListableJobLocator {

    @Autowired
    private ApplicationContext applicationContext;

    public void unregister(String jobName) {

    }

    public void register(JobFactory jobFactory) throws DuplicateJobException {

    }

    public Job getJob(String name) throws NoSuchJobException {
        Job job = (Job) applicationContext.getBean(name);

        if (job == null) {
            throw new NoSuchJobException("IntactContext is not aware of this job: "+name);
        }

        return job;
    }

    public Collection<String> getJobNames() {
        return Arrays.asList(applicationContext.getBeanNamesForType(Job.class));
    }
}
