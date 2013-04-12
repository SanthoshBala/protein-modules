package org.hupo.psi.mi.psicquic.model.server;

import org.apache.http.client.HttpClient;
import org.apache.http.conn.scheme.PlainSocketFactory;
import org.apache.http.conn.scheme.Scheme;
import org.apache.http.conn.scheme.SchemeRegistry;
import org.apache.http.conn.ssl.SSLSocketFactory;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.http.impl.conn.PoolingClientConnectionManager;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.Version;
import org.apache.solr.client.solrj.impl.HttpSolrServer;
import org.mortbay.jetty.Connector;
import org.mortbay.jetty.Server;
import org.mortbay.jetty.nio.SelectChannelConnector;
import org.mortbay.jetty.webapp.WebAppContext;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;

/**
 * solr jetty runner for testing purposes
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>21/06/12</pre>
 */

public class SolrJettyRunner {

    private static Logger log = LoggerFactory.getLogger(SolrJettyRunner.class);

    private int port = 18080;

    protected Server server;
    protected HttpSolrServer solrServer;

    protected File workingDir;
    protected File solrHome;

    public SolrJettyRunner() {
        this(new File(System.getProperty("java.io.tmpdir"), "solr-home-"+System.currentTimeMillis()));
    }

    public SolrJettyRunner(File workingDir) {
        this.workingDir = workingDir;

        if (log.isInfoEnabled()) log.info("Jetty working dir: "+workingDir);

    }

    public void start() throws Exception {
        File solrWar;

        if (workingDir.exists()) {
            solrHome = new File(workingDir, "solr-home");
            solrWar = new File(workingDir, "solr.war");

            if (!solrHome.exists()) {
                throw new IllegalStateException("Working dir "+workingDir+" exists, but no solr-home/ directory could be found");
            }

            if (!solrWar.exists()) {
                throw new IllegalStateException("Working dir "+workingDir+" exists, but no solr.war folder could be found");
            }

            if (log.isDebugEnabled()) log.debug("Using existing directory");

        } else {
            SolrHomeBuilder solrHomeBuilder = new SolrHomeBuilder();

            solrHomeBuilder.install(workingDir);

            solrHome = solrHomeBuilder.getSolrHomeDir();
            solrWar = solrHomeBuilder.getSolrWar();
        }

        // create index folder
        FSDirectory dir = FSDirectory.open(new File(solrHome.getAbsolutePath()+"/data/index"));
        IndexWriterConfig config = new IndexWriterConfig(Version.LUCENE_36, new StandardAnalyzer(Version.LUCENE_36));
        config.setOpenMode(IndexWriterConfig.OpenMode.CREATE);
        IndexWriter writer = new IndexWriter(dir, config);
        writer.commit();
        writer.close();

        System.setProperty("solr.solr.home", solrHome.getAbsolutePath());

        server = new Server();

        Connector connector=new SelectChannelConnector();
        connector.setPort(Integer.getInteger("jetty.port",port).intValue());
        server.setConnectors(new Connector[]{connector});

        WebAppContext webapp = new WebAppContext();
        webapp.setContextPath("/solr");
        webapp.setWar(solrWar.getAbsolutePath());
        webapp.setTempDirectory(workingDir);

        server.setHandler(webapp);

        server.start();
    }

    public void join() throws InterruptedException {
        server.join();
    }

    public void stop() throws Exception {
        if (solrServer != null){
            solrServer.shutdown();
        }
        if (server != null) server.stop();
    }

    public File getSolrHome() {
        return solrHome;
    }

    public String getSolrUrl() {
        return "http://localhost:"+port+"/solr/";
    }

    public HttpSolrServer getSolrServer() {
        solrServer = new HttpSolrServer(getSolrUrl(), createHttpClient());

        solrServer.setConnectionTimeout(5000);
        solrServer.setSoTimeout(5000);
        solrServer.setAllowCompression(true);

        return solrServer;
    }

    protected HttpClient createHttpClient() {
        SchemeRegistry schemeRegistry = new SchemeRegistry();
        schemeRegistry.register(new Scheme("http", 80, PlainSocketFactory
                .getSocketFactory()));
        schemeRegistry.register(new Scheme("https", 443, SSLSocketFactory
                .getSocketFactory()));

        PoolingClientConnectionManager cm = new PoolingClientConnectionManager(schemeRegistry);
        cm.setMaxTotal(128);
        cm.setDefaultMaxPerRoute(24);

        HttpClient httpClient = new DefaultHttpClient(cm);

        return httpClient;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }
}
