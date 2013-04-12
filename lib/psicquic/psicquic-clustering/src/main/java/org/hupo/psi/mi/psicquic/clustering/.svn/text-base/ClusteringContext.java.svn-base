package org.hupo.psi.mi.psicquic.clustering;

import org.hupo.psi.mi.psicquic.clustering.job.dao.ClusteringServiceDaoFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.beans.factory.annotation.Required;
import org.springframework.context.ApplicationContext;
import org.springframework.stereotype.Component;

import java.io.Serializable;

/**
 * Application context.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
@Component
public class ClusteringContext implements Serializable {

    @Autowired
    @Qualifier( value = "inMemory")
    private transient ClusteringServiceDaoFactory daoFactory;

    @Autowired
    private transient ApplicationContext springContext;

    @Autowired
    private ClusteringConfig config;

    public ClusteringContext() {
    }

    ///////////////////////////
    // Getters and Setters

    public ClusteringServiceDaoFactory getDaoFactory() {
        return daoFactory;
    }

    public void setDaoFactory( ClusteringServiceDaoFactory daoFactory ) {
        this.daoFactory = daoFactory;
    }

    public ClusteringConfig getConfig() {
        return config;
    }

    public void setConfig( ClusteringConfig config ) {
        this.config = config;
    }

    public ApplicationContext getSpringContext() {
        return springContext;
    }

    public void setSpringContext( ApplicationContext springContext ) {
        this.springContext = springContext;
    }
}
