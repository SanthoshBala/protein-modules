package org.hupo.psi.mi.psicquic.view.webapp.model;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Controller;

import javax.annotation.PostConstruct;
import javax.faces.context.ExternalContext;
import javax.faces.context.FacesContext;
import javax.faces.event.ActionEvent;
import javax.faces.event.ValueChangeEvent;
import javax.faces.model.SelectItem;
import javax.servlet.http.Cookie;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.Serializable;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 20/08/2012
 * Time: 12:47
 * To change this template use File | Settings | File Templates.
 */
@Controller("columnContext")
@Scope("session")
public class ColumnContextController implements Serializable{

	private static final Log log = LogFactory.getLog(ColumnContextController.class);

	private static String COOKIE_COLS_NAME = "psicquic.cols.view";

	private static String COOKIE_MINIMAL_VALUE = "min cols";
	private static String COOKIE_BASIC_VALUE = "basic cols";
	private static String COOKIE_DEFAULT_VALUE = "def cols";
	private static String COOKIE_ALL_VALUE = "all cols";

	private String[] selectedColumns;
	private List<SelectItem> columnsSelectItems;
	private static String MOLECULE_A_NAME = "moleculeA.name";
	private static String MOLECULE_B_NAME = "moleculeB.name";
	private static String MOLECULE_A_ALTIDS = "moleculeA.altids";
	private static String MOLECULE_B_ALTIDS = "moleculeB.altids";
	private static String MOLECULE_A_ALIASES = "moleculeA.aliases";
	private static String MOLECULE_B_ALIASES = "moleculeB.aliases";
	private static String MOLECULE_A_SPECIES = "moleculeA.species";
	private static String MOLECULE_B_SPECIES = "moleculeB.species";
	private static String FIRST_AUTHOR = "interaction.firstauthor";
	private static String PUBMED_IDENTIFIER = "interaction.pubmedid";
	private static String INTERACTION_TYPE = "interaction.interactiontype";
	private static String INTERACTION_DETECTION_METHOD = "interaction.detectionmethod";
	private static String SOURCE_DATABASE = "interaction.sourcedb";
	private static String INTERACTION_AC = "interaction.ac";
	private static String CONFIDENCE_VALUE = "interaction.confidencevalue";
	private static String EXPANSION_METHOD = "interaction.expansionmethod";
	private static String MOLECULE_A_BIOLOGICAL_ROLE = "moleculeA.biorole";
	private static String MOLECULE_B_BIOLOGICAL_ROLE = "moleculeB.biorole";
	private static String MOLECULE_A_EXPERIMENTAL_ROLE = "moleculeA.exprole";
	private static String MOLECULE_B_EXPERIMENTAL_ROLE = "moleculeB.exprole";
	private static String MOLECULE_A_INTERACTOR_TYPE = "moleculeA.interactortype";
	private static String MOLECULE_B_INTERACTOR_TYPE = "moleculeB.interactortype";
	private static String MOLECULE_A_XREFS = "moleculeA.xrefs";
	private static String MOLECULE_B_XREFS = "moleculeB.xrefs";
	private static String INTERACTION_XREFS = "interaction.xrefs";
	private static String MOLECULE_A_ANNOTATIONS = "moleculeA.annotations";
	private static String MOLECULE_B_ANNOTATIONS = "moleculeB.annotations";
	private static String INTERACTION_ANNOTATIONS = "interaction.annotations";
	private static String HOST_ORGANISM = "interaction.hostorganism";
	private static String INTERACTION_PARAMETERS = "interaction.parameters";
	private static String INTERACTION_CREATION_DATE = "interaction.creationdate";
	private static String INTERACTION_UPDATE_DATE = "interaction.updatedate";
	private static String MOLECULE_A_CHECKSUM = "moleculeA.checksum";
	private static String MOLECULE_B_CHECKSUM = "moleculeB.checksum";
	private static String INTERACTION_CHECKSUM = "interaction.checksum";
	private static String INTERACTION_IS_NEGATIVE = "interaction.negative";
	private static String MOLECULE_A_FEATURES = "moleculeA.features";
	private static String MOLECULE_B_FEATURES = "moleculeB.features";
	private static String MOLECULE_A_STOICHIOMETRY = "moleculeA.stoichiometry";
	private static String MOLECULE_B_STOICHIOMETRY = "moleculeB.stoichiometry";
	private static String MOLECULE_A_PART_IDENT_METHOD = "moleculeA.partidentmethod";
	private static String MOLECULE_B_PART_IDENT_METHOD = "moleculeB.partidentmethod";

	private boolean showTypeRoleIcons;
	private Map<String, Boolean> selectedColumnsMap;

	public ColumnContextController() {
		selectedColumnsMap = Collections.synchronizedMap(new LinkedHashMap<String, Boolean>());
	}

	@PostConstruct
	public void loadColumns() {
		String colsCookie = readCookie(COOKIE_COLS_NAME);

		if (colsCookie != null) {
			if (COOKIE_MINIMAL_VALUE.equals(colsCookie)) {
				selectMinimalColumns();
			} else if (COOKIE_BASIC_VALUE.equals(colsCookie)) {
				selectBasicColumns();
			} else if (COOKIE_DEFAULT_VALUE.equals(colsCookie)) {
				selectDefaultColumns();
			} else if (COOKIE_ALL_VALUE.equals(colsCookie)) {
				selectAllColumns();
			}
		} else {
			selectDefaultColumns();
		}
	}

	private void writeCookie(String name, String value) {
		final ExternalContext externalContext = FacesContext.getCurrentInstance().getExternalContext();

		Cookie cookie = new Cookie(name, value);
		cookie.setMaxAge(3600 * 24 * 200);

		HttpServletResponse response = (HttpServletResponse) externalContext.getResponse();
		response.addCookie(cookie);
	}

	private String readCookie(String name) {
		final ExternalContext externalContext = FacesContext.getCurrentInstance().getExternalContext();
		HttpServletRequest request = (HttpServletRequest) externalContext.getRequest();

		if (request.getCookies() != null) {
			for (Cookie cookie : request.getCookies()) {
				if (name.equals(cookie.getName())) {
					return cookie.getValue();
				}
			}
		}

		return null;
	}

	private String[] getDefaultColumns() {
		return new String[]{
				MOLECULE_A_NAME, MOLECULE_B_NAME, MOLECULE_A_ALIASES, MOLECULE_B_ALIASES, MOLECULE_A_SPECIES,
				MOLECULE_B_SPECIES, INTERACTION_TYPE, FIRST_AUTHOR, PUBMED_IDENTIFIER,CONFIDENCE_VALUE,
				MOLECULE_A_EXPERIMENTAL_ROLE, MOLECULE_B_EXPERIMENTAL_ROLE, INTERACTION_DETECTION_METHOD

		};
	}

	private String[] getBasicColumns() {
		return new String[]{
				MOLECULE_A_NAME, MOLECULE_B_NAME,
				INTERACTION_DETECTION_METHOD, INTERACTION_AC
		};
	}

	private String[] getAllColumns() {
		return new String[]{
				MOLECULE_A_NAME, MOLECULE_B_NAME, MOLECULE_A_ALTIDS, MOLECULE_B_ALTIDS,
				MOLECULE_A_ALIASES, MOLECULE_B_ALIASES, MOLECULE_A_SPECIES, MOLECULE_B_SPECIES, FIRST_AUTHOR, PUBMED_IDENTIFIER,
				INTERACTION_TYPE, INTERACTION_DETECTION_METHOD, SOURCE_DATABASE, INTERACTION_AC, CONFIDENCE_VALUE,
				EXPANSION_METHOD, MOLECULE_A_BIOLOGICAL_ROLE, MOLECULE_B_BIOLOGICAL_ROLE, MOLECULE_A_EXPERIMENTAL_ROLE,
				MOLECULE_B_EXPERIMENTAL_ROLE, MOLECULE_A_INTERACTOR_TYPE, MOLECULE_B_INTERACTOR_TYPE,
				MOLECULE_A_XREFS, MOLECULE_B_XREFS, INTERACTION_XREFS, MOLECULE_A_ANNOTATIONS, MOLECULE_B_ANNOTATIONS,
				INTERACTION_ANNOTATIONS, HOST_ORGANISM, INTERACTION_PARAMETERS, INTERACTION_CREATION_DATE,
				INTERACTION_UPDATE_DATE, MOLECULE_A_CHECKSUM, MOLECULE_B_CHECKSUM, INTERACTION_CHECKSUM,
				INTERACTION_IS_NEGATIVE, MOLECULE_A_FEATURES, MOLECULE_B_FEATURES, MOLECULE_A_STOICHIOMETRY,
				MOLECULE_B_STOICHIOMETRY, MOLECULE_A_PART_IDENT_METHOD, MOLECULE_B_PART_IDENT_METHOD
		};
	}

	private String[] getMinimumColumns() {
		return new String[]{MOLECULE_A_NAME, MOLECULE_B_NAME, INTERACTION_AC};
	}

	public void selectDefaultColumns() {
		selectedColumns = getDefaultColumns();

		selectedColumnsMap.clear();
		for (String columnKey : selectedColumns) {
			selectedColumnsMap.put(columnKey, true);
		}

		writeCookie(COOKIE_COLS_NAME, COOKIE_DEFAULT_VALUE);
	}

	public void selectAllColumns() {
		selectedColumns = getAllColumns();

		selectedColumnsMap.clear();
		for (String columnKey : selectedColumns) {
			selectedColumnsMap.put(columnKey, true);
		}

		writeCookie(COOKIE_COLS_NAME, COOKIE_ALL_VALUE);
	}

	public void selectBasicColumns() {
		this.selectedColumns = getBasicColumns();

		selectedColumnsMap.clear();
		for (String columnKey : selectedColumns) {
			selectedColumnsMap.put(columnKey, true);
		}

		writeCookie(COOKIE_COLS_NAME, COOKIE_BASIC_VALUE);
	}

	public void selectMinimalColumns() {
		this.selectedColumns = getMinimumColumns();
		selectedColumnsMap.clear();
		for (String columnKey : selectedColumns) {
			selectedColumnsMap.put(columnKey, true);
		}

		writeCookie(COOKIE_COLS_NAME, COOKIE_MINIMAL_VALUE);
	}

	private List<SelectItem> createSelectItems() {
		ResourceBundle rb = ResourceBundle.getBundle("org.hupo.psi.mi.psicquic.Messages");

		List<SelectItem> selectItems = new ArrayList<SelectItem>();

		for (String columnKey : getAllColumns()) {
			selectItems.add(new SelectItem(columnKey, rb.getString(columnKey).trim()));
		}

		return selectItems;
	}

	public String[] getSelectedColumns() {
		return selectedColumns;
	}

	public void setSelectedColumns(String[] selectedColumns) {
		if (selectedColumns != null && selectedColumns.length > 0) {
			this.selectedColumns = selectedColumns;
		}
	}

	public List<SelectItem> getColumnsSelectItems() {
		if (columnsSelectItems == null) {
			columnsSelectItems = createSelectItems();
		}
		return columnsSelectItems;
	}

	public void setSelectedColumnsMap(Map<String, Boolean> selectedColumnsMap) {
		this.selectedColumnsMap = selectedColumnsMap;
	}

	public Map<String, Boolean> getSelectedColumnsMap() {
		return selectedColumnsMap;
	}

	public void selectedColumnsChanged(ValueChangeEvent valueChangeEvent) {
		if (valueChangeEvent.getNewValue() instanceof String[]) {
			selectedColumns = ((String[]) valueChangeEvent.getNewValue());
			selectedColumnsMap.clear();
			for (String columnKey : selectedColumns) {
				selectedColumnsMap.put(columnKey, true);
			}
		}
	}

	public void selectedMinimalColumns(ActionEvent actionEvent) {
		selectMinimalColumns();
	}

	public void selectedBasicColumns(ActionEvent actionEvent) {
		selectBasicColumns();
	}

	public void selectedDefaultColumns(ActionEvent actionEvent) {
		selectDefaultColumns();
	}

	public void selectedAllColumns(ActionEvent actionEvent) {
		selectAllColumns();
	}
}
