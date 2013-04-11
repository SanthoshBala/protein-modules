package org.hupo.psi.mi.psicquic.view.webapp.controller.clustering;

import org.hupo.psi.mi.psicquic.model.PsicquicSolrServer;
import org.hupo.psi.mi.psicquic.view.webapp.io.DownloadUtils;
import psidev.psi.mi.tab.model.builder.MitabWriterUtils;
import psidev.psi.mi.tab.model.builder.PsimiTabVersion;

import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 18/09/2012
 * Time: 10:46
 * To change this template use File | Settings | File Templates.
 */
public class ClusteringDownloadUtils extends DownloadUtils {

	private static Map<String, String> contentTypeByFormat;
	private static HashMap<String, String> extensionByFormat;

	private static final String TEXT_PLAIN = "text/plain";
	private static final String TEXT_XML = "text/xml";
	private static final String TXT = "txt";
	private static final String XML = "xml";

	private static final int FILE_NAME_LENGHT = 20;


	public ClusteringDownloadUtils() {
	}

	@Override
	public String getExtension(String format) {
		if (extensionByFormat == null) {
			initializeExtensionByFormat();
		}

		return extensionByFormat.get(format);

	}

	@Override
	public String getContentType(String format) {
		if (contentTypeByFormat == null) {
			initializeContentTypeByFormat();
		}

		return contentTypeByFormat.get(format);
	}

	@Override
	public String getFileName(String query) {
		return "clusteredQuery_" + truncateFileName(query, FILE_NAME_LENGHT) + getDateTime();
	}

	@Override
	public String writeHeaderInMitab(String format) {
		String header = null;

		if(format==null || format.isEmpty()){
			return header;
		}

		if(format.equalsIgnoreCase(PsicquicSolrServer.RETURN_TYPE_MITAB25)){
			header = MitabWriterUtils.buildHeader(PsimiTabVersion.v2_5);
		}
		else if (format.equalsIgnoreCase(PsicquicSolrServer.RETURN_TYPE_MITAB26)){
			header = MitabWriterUtils.buildHeader(PsimiTabVersion.v2_6);
		}
		else if(format.equalsIgnoreCase(PsicquicSolrServer.RETURN_TYPE_MITAB27)){
			header = MitabWriterUtils.buildHeader(PsimiTabVersion.v2_7);
		}

		return header;
	}

	static void initializeContentTypeByFormat() {

		contentTypeByFormat = new HashMap<String, String>();

		contentTypeByFormat.put(PsicquicSolrServer.RETURN_TYPE_XML25, TEXT_XML);
		contentTypeByFormat.put(PsicquicSolrServer.RETURN_TYPE_MITAB25, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSolrServer.RETURN_TYPE_MITAB26, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSolrServer.RETURN_TYPE_MITAB27, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSolrServer.RETURN_TYPE_COUNT, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSolrServer.RETURN_TYPE_DEFAULT, TEXT_PLAIN);

	}

	static void initializeExtensionByFormat() {

		extensionByFormat = new HashMap<String, String>();

		extensionByFormat.put(PsicquicSolrServer.RETURN_TYPE_XML25, XML);
		extensionByFormat.put(PsicquicSolrServer.RETURN_TYPE_MITAB25, TXT);
		extensionByFormat.put(PsicquicSolrServer.RETURN_TYPE_MITAB26, TXT);
		extensionByFormat.put(PsicquicSolrServer.RETURN_TYPE_MITAB27, TXT);
		extensionByFormat.put(PsicquicSolrServer.RETURN_TYPE_COUNT, TXT);
		extensionByFormat.put(PsicquicSolrServer.RETURN_TYPE_DEFAULT, TXT);

	}
}
