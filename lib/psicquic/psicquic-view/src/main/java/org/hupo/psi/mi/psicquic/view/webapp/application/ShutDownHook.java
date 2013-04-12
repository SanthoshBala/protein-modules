package org.hupo.psi.mi.psicquic.view.webapp.application;

import org.apache.myfaces.orchestra.connectionManager.ConnectionManagerDataSource;

import javax.servlet.ServletContextEvent;
import javax.servlet.ServletContextListener;

/**
 * TODO comment this
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>05/12/12</pre>
 */

public class ShutDownHook implements ServletContextListener{
    @Override
    public void contextInitialized(ServletContextEvent sce) {

    }

    @Override
    public void contextDestroyed(ServletContextEvent sce) {

        ConnectionManagerDataSource.releaseAllBorrowedConnections();
        try {
            Thread.sleep(1000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
