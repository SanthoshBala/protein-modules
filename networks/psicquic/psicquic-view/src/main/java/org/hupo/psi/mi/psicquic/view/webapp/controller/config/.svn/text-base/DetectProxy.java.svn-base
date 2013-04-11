package org.hupo.psi.mi.psicquic.view.webapp.controller.config;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.Serializable;
import java.net.InetSocketAddress;
import java.net.Proxy;
import java.net.SocketAddress;

public class DetectProxy implements Serializable {

	Proxy proxy;
	private static final Log log = LogFactory.getLog(DetectProxy.class);

	public DetectProxy() {
		loadProxy();
	}

	private void loadProxy() {

		String proxyHost = System.getProperty("http.proxyHost", null);
		String proxyPort = System.getProperty("http.proxyPort", null);

		if (proxyHost != null) {
			try {
				if (proxyPort != null && proxyPort.length() > 0) {
					int port = Integer.parseInt(proxyPort);
					SocketAddress address = new InetSocketAddress(proxyHost, port);
					setProxy(new Proxy(Proxy.Type.HTTP, address));
				}
			} catch (Exception e) {
				DetectProxy.log.error("Cannot create proxy using port " + proxyPort);
			}
		}
	}

	public Proxy getProxy() {
		return proxy;
	}

	public void setProxy(Proxy proxy) {
		this.proxy = proxy;
	}
}
