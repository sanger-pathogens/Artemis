package uk.ac.sanger.artemis.io;

import java.awt.BorderLayout;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/*
 * A class mediating between the IO package classes and the user. Determines 
 * how the user should be prompted (if at all) during the execution of file 
 * export tasks. It will pop-up windows in swing mode, and will write to 
 * the console in other modes. It will await user interaction in console and
 * swing modes, but will press on in script mode.
 */
public class UI {
	
	/*
	 * Describes 3 possible file export run states.
	 */
	public enum UIMode {
		SWING,		// the traditional GUI application, or running from inside a console with X11 graphics  
		CONSOLE,	// running from a text-only console, able to enter responses into the command-line
		SCRIPT		// running in a completely automated manner, not stopping for errors, useful for cron-jobs
	}
	
	public static UIMode mode = UIMode.SWING;
	
	public static void initalise()
	{
		String uimodeProperty = System.getProperty("uimode");
		if (uimodeProperty != null)
		{
		    uimodeProperty = uimodeProperty.toUpperCase();
			for (UIMode mode : UIMode.values())
			{
				if (mode.toString().equals(uimodeProperty))
				{
					UI.mode = mode;
				}
			}
		}
	}
	
	public static String userInput( final String prompt, final boolean password )
    {
		String result;
		if (password)
		{
			final char[] pChars = System.console().readPassword( prompt + ": " );
			result = new String (pChars);
		} else
		{
			result = System.console().readLine( prompt + ": " );
		}
	    
	    return result.trim();
    }
	
	public static boolean booleanUserInput(final String label, final String message)
	{
		if (mode == UIMode.SCRIPT)
		{
			return true;
		}
		
		else if (mode == UIMode.SWING)
		{
			final JPanel msgPanel = new JPanel(new BorderLayout());
            msgPanel.add(new JLabel(label), BorderLayout.NORTH);
            JTextArea msgError = new JTextArea(message);
            msgError.setLineWrap(true);
            msgError.setEditable(false);
            JScrollPane scollMsg = new JScrollPane(msgError);
            msgPanel.add(scollMsg, BorderLayout.CENTER);
            
            int val = JOptionPane.OK_OPTION;
            if(System.getProperty("noprompt") == null)
              val = JOptionPane.showConfirmDialog(null, msgPanel,
                  "Keys/Qualifier", JOptionPane.OK_CANCEL_OPTION,
                  JOptionPane.QUESTION_MESSAGE);
            
            return (val == JOptionPane.OK_OPTION) ? true : false;
            
		} 
		
		String input = "";
		boolean valid = false;
		while ( valid == false )
		{
			input = userInput(label + "\n" + message + "(y/n)", false);
			if (input.equals("y"))
				valid = true;
			else if (input.equals("n"))
				valid = true;
		}
		if (input.equals("y"))
			return true;
		return false;
	}
	
	public static void info(String message, String heading)
	{
		message("INFO", message, heading);
	}
	
	public static void warn(String message, String heading)
	{
		message("WARN", message, heading);
	}
	
	public static void error(String message, String heading)
	{
		message("ERROR", message, heading);
	}
	
	public static void message(String messageType, String message, String heading)
	{
		switch (mode)
		{
			case SWING:
				JOptionPane.showMessageDialog(null, message,heading, JOptionPane.ERROR_MESSAGE);
			break;
			default:
				System.out.println(messageType + " :: " + heading + " :: " + message);
			break;
		}
	}
	
}
