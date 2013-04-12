package org.hupo.psi.mi.psicquic.registry.config;

import org.apache.axis.utils.XMLUtils;
import org.apache.axis.utils.cache.MethodCache;
import org.apache.cxf.common.util.ReflectionUtil;
import org.springframework.beans.factory.DisposableBean;
import org.springframework.beans.factory.InitializingBean;
import org.springframework.stereotype.Controller;

import java.io.Serializable;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * This class contains some configuration to use multithreading when querying for PSICQUIC services
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>26/03/12</pre>
 */
@Controller
public class PsicquicRegistryThreadConfig implements InitializingBean, DisposableBean, Serializable {

    private ExecutorService executorService;
    private int maxNumberThreads = 50;

    public PsicquicRegistryThreadConfig(){
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
        executorService.shutdown(); // Disable new tasks from being submitted
        try {
            // Wait a while for existing tasks to terminate
            if (!executorService.awaitTermination(30, TimeUnit.SECONDS)) {
                executorService.shutdownNow(); // Cancel currently executing tasks
                // Wait a while for tasks to respond to being cancelled
                if (!executorService.awaitTermination(30, TimeUnit.SECONDS))
                    System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException ie) {
            // (Re-)Cancel if current thread also interrupted
            executorService.shutdownNow();
            // Preserve interrupt status
            Thread.currentThread().interrupt();
        }
    }

    public void afterPropertiesSet() throws Exception {
        executorService = Executors.newFixedThreadPool(maxNumberThreads);
    }

    public void destroy() throws Exception {
        shutDownThreadContext();

        // clear axis thread local
        Method getMc = ReflectionUtil.findMethod(MethodCache.class,
                "getMethodCache");
        Object theObject = MethodCache.getInstance();
        Map mcMap = (Map) getMc.invoke(theObject, null);
        mcMap.clear();

        Field dbfield =
                ReflectionUtil.getDeclaredField(XMLUtils.class,
                        "documentBuilder");
        ThreadLocal db = null;
        try {
            db = (ThreadLocal) dbfield.get(null);
            if (db != null) {
                db.remove();
            }
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
    }
}
