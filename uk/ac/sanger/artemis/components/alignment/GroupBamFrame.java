package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import uk.ac.sanger.artemis.components.Utilities;

class GroupBamFrame extends JFrame
{
  private static final long serialVersionUID = 1L;
  private String[] groups = new String[]{ "Default" };
  private JPanel groupPanel = new JPanel(new GridBagLayout());
  private JPanel bamPanel = new JPanel(new GridBagLayout());
  private HashMap<JCheckBoxMenuItem, JComboBox> bamGroupMap = 
      new HashMap<JCheckBoxMenuItem, JComboBox> ();
  
  private HashMap<String, JCheckBox> groupList = 
      new HashMap<String, JCheckBox>();
  private static final int PAD = 2;
  
  private BamView bamView;
  private JMenu bamFilesMenu;

  GroupBamFrame(final BamView bamView, final JMenu bamFilesMenu)
  {
    super("Alignment Groups");
    
    this.bamView = bamView;
    this.bamFilesMenu = bamFilesMenu;
    
    final JPanel centerPanel = new JPanel(new GridBagLayout());
    final JScrollPane jsp = new JScrollPane(centerPanel);
    getContentPane().add(jsp, BorderLayout.CENTER);
    centerPanel.setBackground(Color.WHITE);

    ///
    /// number of groups
    Box xBox = Box.createHorizontalBox();
    xBox.add(new JLabel("Add Group:"));
    final JTextField newGroup = new JTextField(10);
    xBox.add(newGroup);
    newGroup.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        addGroup(bamView, bamFilesMenu, newGroup, jsp);
      }
    });
    newGroup.addFocusListener(new java.awt.event.FocusAdapter() {
      public void focusGained(java.awt.event.FocusEvent evt) {
        SwingUtilities.invokeLater( new Runnable() {
          public void run() {
            newGroup.selectAll();              
          }
        });
      }
    });

    final JButton addGroupButton = new JButton("ADD");
    xBox.add(addGroupButton);
    addGroupButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        if(newGroup.getText() == null ||
           newGroup.getText().trim().equals(""))
        {
          JOptionPane.showMessageDialog(GroupBamFrame.this, 
              "Type in the name for the new group.", "Group name missing", 
              JOptionPane.INFORMATION_MESSAGE);
          return;
        }
        addGroup(bamView, bamFilesMenu, newGroup, jsp);
      }
    });
    xBox.add(Box.createHorizontalGlue());
    getContentPane().add(xBox, BorderLayout.NORTH);
    
    //
    xBox = Box.createHorizontalBox();
    final JButton close = new JButton("CLOSE");
    close.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        setVisible(false);
      }
    });
    xBox.add(close);
    xBox.add(Box.createHorizontalGlue());
    getContentPane().add(xBox, BorderLayout.SOUTH); 

    
    ///
    ///
    final GridBagConstraints c = new GridBagConstraints();
    c.gridy = 0;
    c.gridx = 0;
    c.ipadx = PAD;
    c.ipady = PAD;
    c.anchor = GridBagConstraints.WEST;
    c.gridwidth = 2;
    final JLabel alignLabel = new JLabel("ALIGNMENTS:");
    alignLabel.setFont(alignLabel.getFont().deriveFont(Font.BOLD));
    centerPanel.add(alignLabel, c);
    c.gridwidth = 1;
    c.gridy += 1;
    centerPanel.add(Box.createVerticalStrut(5), c);
    c.gridy += 1;

    bamPanel.setBackground(Color.WHITE);
    centerPanel.add(bamPanel, c);
    
    c.gridx = 10;
    c.gridheight = GridBagConstraints.REMAINDER;
    centerPanel.add(groupPanel, c);
    groupPanel.setBackground(Color.WHITE);

    Utilities.rightJustifyFrame(this);
  }
  
  protected void updateAndDisplay()
  {
    updateBamPanel();
    updateGroupPanel();
    pack();
    setVisible(true);
  }
  
  private void addGroup(final BamView bamView, final JMenu bamFilesMenu,
                        final JTextField newGroup, final JScrollPane jsp)
  {
    final String tmpGroups[] = new String[groups.length+1];
    for(int i=0; i<groups.length; i++)
      tmpGroups[i] = groups[i];
    tmpGroups[tmpGroups.length-1] = newGroup.getText().trim();
    groups = tmpGroups;
    
    updateGroupPanel();
    updateBamPanel();
    jsp.validate();
  }
  
  private void updateBamPanel()
  {
    final Component cs[] = bamFilesMenu.getMenuComponents();
    bamPanel.removeAll();
    final GridBagConstraints c = new GridBagConstraints();
    c.gridy = 0;
    c.gridx = 0;
    c.ipadx = PAD;
    c.ipady = PAD;
    c.anchor = GridBagConstraints.WEST;
    
    for(Component cp : cs)
    {
      if(cp instanceof JCheckBoxMenuItem)
      {
        final JCheckBoxMenuItem cbBam = (JCheckBoxMenuItem) cp;
        final String bam = cbBam.getText();
        c.gridy += 1;
        c.gridx = 0;
        bamPanel.add(new JLabel(bam), c);
        
        Color col = bamView.getColorByJCheckBoxMenuItem(cbBam);
        if(col == null)
          col = Color.BLACK;
        
        c.gridx = 1;
        bamPanel.add(new JLabel(bamView.getImageIcon(col)), c);

        c.gridx = 2;
        JComboBox groupList = new JComboBox( groups );
        bamPanel.add(groupList, c);
        
        if(bamGroupMap.containsKey(cbBam))
          groupList.setSelectedItem(
              (String) bamGroupMap.get(cbBam).getSelectedItem());

        bamGroupMap.put(cbBam, groupList);
      }
    }
  }
  
  private void updateGroupPanel()
  {
    groupPanel.removeAll();
    final GridBagConstraints c = new GridBagConstraints();
    c.gridy = 0;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    
    final JLabel groupLabel = new JLabel("GROUPS:");
    groupLabel.setFont(groupLabel.getFont().deriveFont(Font.BOLD));
    groupPanel.add(groupLabel, c);
    c.gridy += 2;
    for(final String group: groups)
    {
      c.gridy += 1;
      
      final JCheckBox cbox;
      if(groupList.containsKey(group))
        cbox = groupList.get(group);
      else
      {
        cbox = new JCheckBox(group, true);
        groupList.put(group, cbox);
      }
      groupPanel.add(cbox, c);
      cbox.addItemListener(new ItemListener(){
        public void itemStateChanged(ItemEvent arg0)
        {
          Iterator<JCheckBoxMenuItem> it = bamGroupMap.keySet().iterator();
          while(it.hasNext())
          {
            JCheckBoxMenuItem cbs = it.next();
            String thisGroup = (String) bamGroupMap.get(cbs).getSelectedItem();
            if(group.equals(thisGroup))
              cbs.setSelected(cbox.isSelected());
          }
        }
      });
    }
  }
  

}