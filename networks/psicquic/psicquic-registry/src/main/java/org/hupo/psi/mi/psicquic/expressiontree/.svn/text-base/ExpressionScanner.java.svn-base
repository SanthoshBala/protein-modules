package org.hupo.psi.mi.psicquic.expressiontree;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * 
 */

/**
 * @author Erik Pfeiffenberger
 *
 */
public class ExpressionScanner {
	private String [] characters;
	private String text;
	
	private int currentCharPosition;
	
	private int currentSymbolPosition;
	private List<String> symbolsList;
	
	private final String NOT_SYMBOL = "not";
	private final String AND_SYMBOL = "and";
	private final String OR_SYMBOL = "or";
	private final String LEFT_PAR_SYMBOL = "(";
	private final String RIGHT_PAR_SYMBOL = ")";
	
	private ExpressionStates currState = ExpressionStates.FreeText;
	
	private int charsRead = 0;
	
	private final Pattern NOT_SYMBOL_PATTERN = Pattern.compile("^"+NOT_SYMBOL);
	private final Pattern AND_SYMBOL_PATTERN = Pattern.compile("^"+AND_SYMBOL);
	private final Pattern OR_SYMBOL_PATTERN = Pattern.compile("^"+OR_SYMBOL);
	
	private final Pattern LEFT_PAR_SYMBOL_PATTERN = Pattern.compile("^\\"+LEFT_PAR_SYMBOL);
	private final Pattern RIGHT_PAR_SYMBOL_PATTERN = Pattern.compile("^\\"+RIGHT_PAR_SYMBOL);
	
	private final Pattern MI_ID_SYMBOL_PATTERN = Pattern.compile("^mi:\\d{4}");
	
	private final Pattern FREE_TEXT_PATTERN = Pattern.compile("^[a-zA-Z0-9\\s]+");
	

	
	private enum ExpressionStates{
		FreeText
	}
	
	
	public ExpressionScanner(String text) {
		
		text = text.toLowerCase();
		this.text = text;
		String [] tmp = text.split("");

		characters = new String [tmp.length -1];
		
		//to remove "" on position 0
		for (int i=1; i < tmp.length; i++){
			characters[i-1] = tmp[i];
		}
		
		currentSymbolPosition = -1;
		symbolsList = new ArrayList<String>();
		
		currentCharPosition = 0;
	}
	
	public String nextSymbol(){
		text = text.trim();
		
		String symbol = "";
		
		if(text.startsWith(RIGHT_PAR_SYMBOL)){
			symbol = RIGHT_PAR_SYMBOL;
			text = text.substring(RIGHT_PAR_SYMBOL.length());
		} else if (text.startsWith(LEFT_PAR_SYMBOL)){
			symbol = LEFT_PAR_SYMBOL;
			text = text.substring(LEFT_PAR_SYMBOL.length());
		}else if (text.startsWith(AND_SYMBOL+" ") || text.startsWith(AND_SYMBOL+"(")){
			symbol = AND_SYMBOL;
			text = text.substring(AND_SYMBOL.length());
		}else if (text.startsWith(OR_SYMBOL+" ") || text.startsWith(OR_SYMBOL+"(")){
			symbol = OR_SYMBOL;
			text = text.substring(OR_SYMBOL.length());
		}else if (text.startsWith(NOT_SYMBOL+" ") || text.startsWith(NOT_SYMBOL+"(")){
			symbol = NOT_SYMBOL;
			text = text.substring(NOT_SYMBOL.length());
		} else if(text.matches("^\".+\"")){
			text = text.substring(1);
			int index = text.indexOf("\"");
			symbol = text.substring(0, index);
			text = text.substring(index+1);
		}else if(text.matches("^mi:\\d{4}")){
			
			//every mi identifier is 7 characters long
			symbol = text.substring(0, 7);
			//every mi identifier is 7 characters long
			text = text.substring(7);
		} else {
			String [] arr = text.split("( and )|( and\\()|( and\\))|( or )|( or\\()|( or\\))|( not )|( not\\()|( not\\))|(\\()|(\\))");
			String firstElement = arr[0];
			symbol = firstElement;
			text = text.substring(firstElement.length());
			
		}
		
		symbol = symbol.trim();
		currentSymbolPosition++;
		symbolsList.add(symbol);
		
		return symbol;
	}
	
	public boolean isNot(){
		return currentSymbolPosition == -1? false : NOT_SYMBOL_PATTERN.matcher(symbolsList.get(currentSymbolPosition)).matches();
		
	}
	
	public boolean isAnd(){
		return currentSymbolPosition == -1? false : AND_SYMBOL_PATTERN.matcher(symbolsList.get(currentSymbolPosition)).matches();
		
	}
	
	public boolean isOr(){
		return currentSymbolPosition == -1? false : OR_SYMBOL_PATTERN.matcher(symbolsList.get(currentSymbolPosition)).matches();
		
	}
	
	public boolean isLeftPar(){
		return currentSymbolPosition == -1? false : LEFT_PAR_SYMBOL_PATTERN.matcher(symbolsList.get(currentSymbolPosition)).matches();
		
	}
	
	public boolean isRightPar(){
		return currentSymbolPosition == -1? false : RIGHT_PAR_SYMBOL_PATTERN.matcher(symbolsList.get(currentSymbolPosition)).matches();
		
	}
	
	public boolean isMiIdentifier(){
		return currentSymbolPosition == -1? false : MI_ID_SYMBOL_PATTERN.matcher(symbolsList.get(currentSymbolPosition)).matches();
		
	}
	
	public boolean isFreeText(){
		boolean isFreeText = false;
		
		if(currentSymbolPosition == -1){
			isFreeText = false;
		} else if(!isAnd() && !isOr() && !isNot() && !isMiIdentifier() && !isLeftPar() && !isRightPar()){
			
			String currSymbol = symbolsList.get(currentSymbolPosition);
			
			if (currSymbol.matches("^.+")){
				isFreeText = true;
			} else {
				isFreeText = false;
			}
			
			
		}
		
		return isFreeText;
		
		
	}
	
	public boolean hasNext(){
		return text.length() > 0;
	}
	
	public String getCurrSymbol(){
		String currSymbol = "";
		
		if(currentSymbolPosition >= 0){
			currSymbol = symbolsList.get(currentSymbolPosition);
			
		}
		
		return currSymbol;
	}
	

}
