package org.hupo.psi.mi.psicquic.model.server;

import org.apache.commons.io.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Build solr home
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>21/06/12</pre>
 */

public class SolrHomeBuilder {

    private static Logger log = LoggerFactory.getLogger(SolrHomeBuilder.class);

    private File solrHomeDir;
    private File solrWar;

    public SolrHomeBuilder() {
    }

    public void install(File solrWorkingDir) throws IOException {
        if (log.isInfoEnabled()) log.info("Installing Intact SOLR Home at: "+solrWorkingDir);

        // copy resource directory containing solr-home and war file
        File solrHomeToCreate = new File(solrWorkingDir,"solr-home");
        File solrWarToCreate = new File(solrWorkingDir,"solr.war");

        // only copy solr-home when solr-home does not exist
        if (!solrHomeToCreate.exists() && solrHomeDir == null){
            solrHomeToCreate.mkdirs();

            File solrHomeToCopy = new File(SolrHomeBuilder.class.getResource("/solr-home").getFile());
            // is in the resources
            if (solrHomeToCopy.exists()){
                FileUtils.copyDirectory(solrHomeToCopy, solrHomeToCreate);

                if (!solrWarToCreate.exists() && solrWar == null){
                    File solrWarToCopy = new File(SolrHomeBuilder.class.getResource("/solr.war").getFile());
                    FileUtils.copyFile(solrWarToCopy, solrWarToCreate);
                }
            }
            // is in the jar in the dependencies
            else {

                String originalName = SolrHomeBuilder.class.getResource("/solr-home").getFile();
                String jarFileName = originalName.substring(0, originalName.indexOf("!")).replace("file:", "");

                JarFile jarFile = new JarFile(jarFileName);
                Enumeration<JarEntry> jarEntries = jarFile.entries();

                // write
                while (jarEntries.hasMoreElements()) {
                    JarEntry entry = jarEntries.nextElement();

                    // solr war file
                    if (entry.getName().endsWith("solr.war") && !solrWarToCreate.exists() && solrWar == null) {

                        InputStream inputStream = jarFile.getInputStream(entry);

                        try{
                            FileUtils.copyInputStreamToFile(inputStream, solrWarToCreate);
                        }
                        finally {
                            inputStream.close();
                        }
                    }
                    else if (entry.toString().startsWith("solr-home")){
                        File fileToCreate = new File(solrWorkingDir, entry.toString());

                        if (entry.isDirectory()) {
                            fileToCreate.mkdirs();
                            continue;
                        }

                        InputStream inputStream = jarFile.getInputStream(entry);

                        try{
                            FileUtils.copyInputStreamToFile(inputStream, fileToCreate);
                        }
                        finally {
                            inputStream.close();
                        }
                    }
                }
            }

            solrHomeDir = solrHomeToCreate;
            solrWar = solrWarToCreate;
        }
        // only copy solr.war when solr.war does not exist
        else if (!solrWarToCreate.exists() && solrWar == null){

            File solrHomeToCopy = new File(SolrHomeBuilder.class.getResource("/solr-home").getFile());
            // is in the resources
            if (solrHomeToCopy.exists()){
                File solrWarToCopy = new File(SolrHomeBuilder.class.getResource("/solr.war").getFile());
                FileUtils.copyFileToDirectory(solrWarToCopy, solrWorkingDir);
            }
            // is in the jar in the dependencies
            else {

                String originalName = SolrHomeBuilder.class.getResource("/solr-home").getFile();
                String jarFileName = originalName.substring(0, originalName.indexOf("!")).replace("file:", "");

                JarFile jarFile = new JarFile(jarFileName);
                Enumeration<JarEntry> jarEntries = jarFile.entries();

                // write
                while (jarEntries.hasMoreElements()) {
                    JarEntry entry = jarEntries.nextElement();

                    // solr war file
                    if (entry.getName().endsWith("solr.war")) {
                        File fileToCreate = new File(solrWorkingDir, entry.toString());

                        InputStream inputStream = jarFile.getInputStream(entry);

                        try{
                            FileUtils.copyInputStreamToFile(inputStream, fileToCreate);
                        }
                        finally {
                            inputStream.close();
                        }
                    }
                }
            }

            solrHomeDir = solrHomeToCreate;
            solrWar = solrWarToCreate;
        }

        if (log.isDebugEnabled()) {
            log.debug("\nSolr Home: {}\nSolr WAR: {}", solrHomeDir.toString(), solrWar.toString());
        }

    }

    public File getSolrHomeDir() {
        return solrHomeDir;
    }

    public File getSolrWar() {
        return solrWar;
    }

    public void setSolrHomeDir(File solrHomeDir) {
        this.solrHomeDir = solrHomeDir;
    }

    public void setSolrWar(File solrWar) {
        this.solrWar = solrWar;
    }
}
