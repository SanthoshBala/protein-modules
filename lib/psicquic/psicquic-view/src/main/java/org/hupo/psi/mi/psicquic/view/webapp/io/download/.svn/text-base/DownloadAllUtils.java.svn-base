package org.hupo.psi.mi.psicquic.view.webapp.io.download;

import org.hupo.psi.mi.psicquic.view.webapp.io.DownloadUtils;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;
import psidev.psi.mi.tab.model.builder.MitabWriterUtils;
import psidev.psi.mi.tab.model.builder.PsimiTabVersion;

import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 30/10/2012
 * Time: 16:15
 * To change this template use File | Settings | File Templates.
 */
public class DownloadAllUtils extends DownloadUtils {

	private static Map<String, String> contentTypeByFormat;
	private static HashMap<String, String> extensionByFormat;
	private static final String TEXT_PLAIN = "text/plain";
	private static final String TXT = "txt";

	private static final int FILE_NAME_LENGHT = 20;

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
		return "results_" + truncateFileName(query, FILE_NAME_LENGHT) + getDateTime();
	}

	@Override
	public String writeHeaderInMitab(String format) {
		String header = null;

		if(format==null || format.isEmpty()){
			return header;
		}

		if(format.equalsIgnoreCase(PsicquicSimpleClient.MITAB25)){
			header = MitabWriterUtils.buildHeader(PsimiTabVersion.v2_5);
		}
		else if (format.equalsIgnoreCase(PsicquicSimpleClient.MITAB26)){
			header = MitabWriterUtils.buildHeader(PsimiTabVersion.v2_6);
		}
		else if(format.equalsIgnoreCase(PsicquicSimpleClient.MITAB27)){
			header = MitabWriterUtils.buildHeader(PsimiTabVersion.v2_7);
		}

		return header;
	}

	static void initializeContentTypeByFormat() {

		contentTypeByFormat = new HashMap<String, String>();

		contentTypeByFormat.put(PsicquicSimpleClient.MITAB25, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSimpleClient.MITAB26, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSimpleClient.MITAB27, TEXT_PLAIN);
		contentTypeByFormat.put(PsicquicSimpleClient.COUNT, TEXT_PLAIN);

	}

	static void initializeExtensionByFormat() {

		extensionByFormat = new HashMap<String, String>();

		extensionByFormat.put(PsicquicSimpleClient.MITAB25, TXT);
		extensionByFormat.put(PsicquicSimpleClient.MITAB26, TXT);
		extensionByFormat.put(PsicquicSimpleClient.MITAB27, TXT);
		extensionByFormat.put(PsicquicSimpleClient.COUNT, TXT);

	}
}
