package org.hupo.psi.mi.psicquic.view.webapp.application;

import org.springframework.beans.factory.DisposableBean;
import org.springframework.beans.factory.InitializingBean;
import org.springframework.stereotype.Controller;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * This class contains some configuration to use multithreading when querying for PSICQUIC services
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>26/03/12</pre>
 */
@Controller
public class PsicquicThreadConfig implements InitializingBean, DisposableBean {

    private ExecutorService executorService;
    private int maxNumberThreads = 250;

    public PsicquicThreadConfig(){
    }

    public ExecutorService getExecutorService() {
        return executorService;
    }

    public void setExecutorService(ExecutorService executorService) {
        this.executorService = executorService;
    }

    public int getMaxNumberThreads() {
        return maxNumberThreads;
    }

    public void setMaxNumberThreads(int maxNumberThreads) {
        this.maxNumberThreads = maxNumberThreads;
    }

    public void shutDownThreadContext(){
        executorService.shutdownNow();
    }

    @Override
    public void afterPropertiesSet() throws Exception {
        executorService = Executors.newFixedThreadPool(maxNumberThreads);
    }

    @Override
    public void destroy() throws Exception {
        shutDownThreadContext();
    }
}
