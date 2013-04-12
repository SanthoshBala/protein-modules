/*
 * Copyright (c) 2002 The European Bioinformatics Institute, and others.
 * All rights reserved. Please see the file LICENSE
 * in the root directory of this distribution.
 */
package org.hupo.psi.mi.psicquic.view.webapp;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.springframework.beans.factory.InitializingBean;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationContext;
import org.springframework.stereotype.Component;

import java.util.EventListener;

/**
 * Initializes the application.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
@Component
public class ApplicationInitializer implements InitializingBean, EventListener {

    public static final Log log = LogFactory.getLog( ApplicationInitializer.class);

    @Autowired
    private ApplicationContext springContext;

    @Autowired
    private ClusteringContext clusteringContext;

    public ApplicationInitializer() {
    }

    public void afterPropertiesSet() throws Exception {

        if( springContext != null ) {
            System.out.println("Location (clusteringContext.getConfig().getDataLocation()): " + clusteringContext.getConfig().getDataLocation());

            if( clusteringContext.getSpringContext() == null ) {
                log.error( "Failed to initialize ClusteringContext: springContext is null." );
            }
        } else {
            log.error( "||||||||||||||||||||||||||| SPRING WAS NOT INITIALIZED ||||||||||||||||||||||" );
        }
    }
}
