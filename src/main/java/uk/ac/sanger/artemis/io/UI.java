package uk.ac.sanger.artemis.io;

import java.awt.BorderLayout;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.IOException;

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
	
	private static String input()
	{
	    InputStreamReader reader = new InputStreamReader (System.in);
        BufferedReader buf_reader = new BufferedReader (reader);
        String result = "";
        try {
          result = buf_reader.readLine().trim(); // read the number as a string
        }
        catch (IOException ioe) {
          System.out.println ("IO exception = " + ioe);
        }
        return result.trim();
	}
	
	public static String userInput( final String prompt, final boolean password )
    {
        if (password == true)
        {
            UIEraserThread et = new UIEraserThread(prompt + ": ");
            Thread mask = new Thread(et);
            mask.start();
            String result = input();
            et.stopMasking();
            mask.stop();
            // System.out.println("Password: '" + result + "'");
            return result;
        }
        System.out.print(prompt + ": " );
        return input();
    }
	
	public static boolean booleanUserInput(final String label, final String message)
	{
		if (mode == UIMode.SCRIPT)
		{
		    System.out.println(label + " : " + message + " (in script mode so continuing...)");
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
            //if(System.getProperty("noprompt") == null)
            val = JOptionPane.showConfirmDialog(null, msgPanel,
                  "Keys/Qualifier", JOptionPane.OK_CANCEL_OPTION,
                  JOptionPane.QUESTION_MESSAGE);
            
            return (val == JOptionPane.OK_OPTION) ? true : false;
            
		} 
		
		return boolConsoleInput(label + "\n" + message + "(y/n)");
	}
	
	private static boolean boolConsoleInput(String prompt)
	{
	    System.out.print(prompt + ": " );
	    InputStreamReader reader = new InputStreamReader (System.in);
        BufferedReader buf_reader = new BufferedReader (reader);
        String result = "";
        try {
          result = buf_reader.readLine().trim(); // read the number as a string
        }
        catch (IOException ioe) {
          System.out.println ("IO exception = " + ioe);
        }
	    // System.out.println("Entered '" + result + "'");
	    if (result.equals("y"))
	        return true;
	    if (result.equals("n"))
	        return false;
	    System.out.println("Please answer y or n");
	    return boolConsoleInput(prompt);
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
	    if ((mode == UIMode.SCRIPT) || (mode == UIMode.CONSOLE))
        {
            System.out.println(messageType + " :: " + heading + " :: " + message);
	    }
	    else if (mode == UIMode.SWING)
	    {
	        JOptionPane.showMessageDialog(null, message,heading, JOptionPane.ERROR_MESSAGE);
	    } 
	}
	
}
