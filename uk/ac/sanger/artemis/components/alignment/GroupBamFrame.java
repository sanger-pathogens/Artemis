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
        addGroup(newGroup.getText().trim());
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
        addGroup(newGroup.getText().trim());
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
  
  protected void addGroup(final String newGroup)
  {
    final String tmpGroups[] = new String[groups.length+1];
    System.arraycopy(groups, 0, tmpGroups, 0, groups.length);
    tmpGroups[tmpGroups.length-1] = newGroup;
    groups = tmpGroups;
    
    updateGroupPanel();
    updateBamPanel();
    validate();
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

        c.gridx = 1;
        bamPanel.add(new JLabel(cbBam.getIcon()), c);

        c.gridx = 2;
        JComboBox groupCombo = new JComboBox( groups );
        bamPanel.add(groupCombo, c);
        
        if(bamGroupMap.containsKey(cbBam))
          groupCombo.setSelectedItem(
              (String) bamGroupMap.get(cbBam).getSelectedItem());

        bamGroupMap.put(cbBam, groupCombo);
      }
    }
  }
  
  /**
   * Add a BAM to a specified group
   * @param bam
   * @param group
   */
  protected void addToGroup(String bam, String group)
  {
    final Component cs[] = bamFilesMenu.getMenuComponents();
    for(Component cp : cs)
    {
      if(cp instanceof JCheckBoxMenuItem)
      {
        final JCheckBoxMenuItem cbBam = (JCheckBoxMenuItem) cp;
        final String thisBam = cbBam.getText();
        if(bam.equals(thisBam))
        {
          bamGroupMap.get(cbBam).setSelectedItem(group);
          return;
        }
      }
    }
  }
  
  /**
   * For a give BAM file return the group it belongs to.
   * @param bamName
   * @return
   */
  protected String getGroupName(final String bamName)
  {
    Iterator<JCheckBoxMenuItem> it = bamGroupMap.keySet().iterator();
    while(it.hasNext())
    {
      JCheckBoxMenuItem cbs = it.next();
      if(cbs.getText().equals(bamName))
        return (String) bamGroupMap.get(cbs).getSelectedItem();
    }
    return null;
  }
  
  protected int getNumberOfGroups()
  {
    return groups.length;
  }
  
  /**
   * Find the maximum number of BAM files found in the groups.
   * @return
   */
  protected int getMaximumBamsInGroup()
  {
    int max = 1;
    int grpMax[] = new int[getNumberOfGroups()];
    for(int i=0; i<grpMax.length; i++)
      grpMax[i] = 0;
    Iterator<JCheckBoxMenuItem> it = bamGroupMap.keySet().iterator();
    while(it.hasNext())
    {
      JCheckBoxMenuItem cbs = it.next();
      String grp = (String) bamGroupMap.get(cbs).getSelectedItem();
      for(int i=0; i<groups.length; i++)
      {
        if(grp.equals(groups[i]))
        {
          grpMax[i]++;
          break;
        }
      }
    }
    
    for(int i=0; i<grpMax.length; i++)
    {
      if(grpMax[i] > max)
        max = grpMax[i];
    }
    return max;
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