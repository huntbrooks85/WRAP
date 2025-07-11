// ----------------------------------------------------------------------- //
//  WRAP.java v2.1.0
//  By Hunter Brooks, at NAU/UToledo, Flagstaff: April 23, 2025
// ----------------------------------------------------------------------- //



// Import All Packages
// ------------------------------------------------------------- //
// Import AWT Packages
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
// ------------------------------------------------------------- //



// Makes Main GUI Class
// ------------------------------------------------------------- //
public class WRAP {
    public static void main(String[] args) {
        // Frame Setup
        // ------------------------------------------------------------- //
        JFrame frame = new JFrame("WRAP");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(520, 255);
        frame.setLayout(null);
        frame.setResizable(false);

        // Finds Path of File
        String jarPath = "";
        try {
            jarPath = new File(WRAP.class.getProtectionDomain().getCodeSource().getLocation().toURI()).getPath();
        } catch (URISyntaxException e) {
            e.printStackTrace();
        }
        File appDir = new File(jarPath).getParentFile();
        // ------------------------------------------------------------- //



        // File Loader Panel
        // ------------------------------------------------------------- //
        FileLoaderPanel fileLoaderPanel = new FileLoaderPanel();
        fileLoaderPanel.setBounds(10, 10, 150, 220); 
        frame.add(fileLoaderPanel);
        // ------------------------------------------------------------- //



        // WRAP Data Panel
        // ------------------------------------------------------------- //
        JPanel mainControlPanel = new JPanel();
        mainControlPanel.setLayout(null);
        mainControlPanel.setPreferredSize(new Dimension(200, 350));
        mainControlPanel.setBounds(165, 10, 345, 200);
        mainControlPanel.setBorder(BorderFactory.createTitledBorder("Input Requirements"));

        // RA, DEC, RADIUS Input Panel
        // JPanel coordsPanel = new JPanel(new GridLayout(1, 3));

        // Creates RA Field with grey text in box
        JTextField raField = new JTextField("R.A. (deg)", 20);
        raField.setForeground(Color.GRAY);
        raField.setBounds(10, 20, 80, 20);
        raField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (raField.getText().equals("R.A. (deg)")) {
                    raField.setText("");
                    raField.setForeground(Color.BLACK);
                }
            }

            @Override
            public void focusLost(FocusEvent e) {
                if (raField.getText().isEmpty()) {
                    raField.setText("R.A. (deg)");
                    raField.setForeground(Color.GRAY);
                }
            }
        });

        // Creates DEC Field with grey text in box
        JTextField decField = new JTextField("Decl. (deg)", 20);
        decField.setBounds(90, 20, 80, 20);
        decField.setForeground(Color.GRAY);
        decField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (decField.getText().equals("Decl. (deg)")) {
                    decField.setText("");
                    decField.setForeground(Color.BLACK);
                }
            }

            @Override
            public void focusLost(FocusEvent e) {
                if (decField.getText().isEmpty()) {
                    decField.setText("Decl. (deg)");
                    decField.setForeground(Color.GRAY);
                }
            }
        });

        // Creates Radius Field with grey text in box
        JTextField radiusField = new JTextField("Search Radius (arcsec)", 20);
        radiusField.setForeground(Color.GRAY);
        radiusField.setBounds(170, 20, 170, 20);
        radiusField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (radiusField.getText().equals("Search Radius (arcsec)")) {
                    radiusField.setText("");
                    radiusField.setForeground(Color.BLACK);
                }
            }

            @Override
            public void focusLost(FocusEvent e) {
                if (radiusField.getText().isEmpty()) {
                    radiusField.setText("Search Radius (arcsec)");
                    radiusField.setForeground(Color.GRAY);
                }
            }
        });

        // Creates a "Multiquery" checkbox area
        JCheckBox checkboxmultiquery = new JCheckBox("Multi-Query");
        checkboxmultiquery.setFont(new Font("Times New Roman", Font.PLAIN, 9));
        checkboxmultiquery.setBounds(2, 155, 80, 20);
        mainControlPanel.add(checkboxmultiquery);

        // Listener to swap components
        JTextField fileField = new JTextField("path/to/file.csv");
        fileField.setForeground(Color.GRAY);
        fileField.setBounds(10, 20, 150, 20);
        fileField.setEditable(false);
        fileField.setVisible(false);
        
        JButton browseButton = new JButton("Browse");
        browseButton.setBounds(160, 20, 50, 20);
        browseButton.setVisible(false);
        
        browseButton.addActionListener(e -> {
            JFileChooser fileChooser = new JFileChooser();
            int result = fileChooser.showOpenDialog(mainControlPanel);
            if (result == JFileChooser.APPROVE_OPTION) {
                if (fileField.getText().toLowerCase().endsWith(".csv")) {
                    fileField.setText(fileChooser.getSelectedFile().getAbsolutePath());
                }else{
                    JOptionPane.showMessageDialog(fileField, "Input a CSV File", "Input Error", JOptionPane.ERROR_MESSAGE);
                }
            }
        });
        
        checkboxmultiquery.addActionListener(e -> {
            boolean selected = checkboxmultiquery.isSelected();
            raField.setVisible(!selected);
            decField.setVisible(!selected);
            fileField.setVisible(selected);
            browseButton.setVisible(selected);
        
            mainControlPanel.revalidate();
            mainControlPanel.repaint();
        });        

        // Adds RA, DECL, RADIUS, and FILE Fields
        mainControlPanel.add(raField);
        mainControlPanel.add(decField);
        mainControlPanel.add(radiusField);
        mainControlPanel.add(fileField);
        mainControlPanel.add(browseButton);

        // Catalog Panel with checkboxes
        JPanel catalogPanel = new JPanel(new GridLayout(0, 2));
        JCheckBox checkBox1 = new JCheckBox("CatWISE");
        JCheckBox checkBox2 = new JCheckBox("AllWISE");
        JCheckBox checkBox3 = new JCheckBox("Gaia");
        JCheckBox checkBox4 = new JCheckBox("VHS");
        JCheckBox checkBox5 = new JCheckBox("UKIDSS");
        JCheckBox checkBox6 = new JCheckBox("2MASS");
        JCheckBox checkBox7 = new JCheckBox("SDSS");
        JCheckBox checkBox8 = new JCheckBox("PanSTARRS");
        JCheckBox checkBox9 = new JCheckBox("NSC");
        JCheckBox checkBox10 = new JCheckBox("GALEX");

        // Set horizontal alignment for the text in the checkboxes
        JCheckBox[] checkboxes = {checkBox1, checkBox2, checkBox3, checkBox4, checkBox5, checkBox6, checkBox7, checkBox8, checkBox9, checkBox10};
        for (JCheckBox checkBox : checkboxes) {
            checkBox.setFont(new Font("Times New Roman", Font.PLAIN, 16));
            catalogPanel.add(checkBox);
        }

        // Add checkboxes to the catalog panel
        JScrollPane scrollPanelcatalog = new JScrollPane(catalogPanel);
        scrollPanelcatalog.setBounds(8, 45, 330, 100);
        mainControlPanel.add(scrollPanelcatalog);

        // Select/Deselect/Other
        JPanel optionsPanel = new JPanel();
        JButton checkboxall = new JButton("Select All");
        checkboxall.setPreferredSize(new Dimension(70, 15));
        checkboxall.setFont(new Font("Times New Roman", Font.PLAIN, 10));
        JButton checkboxdeall = new JButton("Deselect All");
        checkboxdeall.setPreferredSize(new Dimension(70, 15));
        checkboxdeall.setFont(new Font("Times New Roman", Font.PLAIN, 10));
        optionsPanel.add(checkboxall);
        optionsPanel.add(checkboxdeall);
        optionsPanel.setBounds(100, 145, 150, 20);
        mainControlPanel.add(optionsPanel);

        // Array or list to hold all catalog checkboxes
        JCheckBox[] catalogCheckboxes = {
            checkBox1, checkBox2, checkBox3, checkBox4, checkBox5,
            checkBox6, checkBox7, checkBox8, checkBox9, checkBox10
        };

        // Select All listener
        checkboxall.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (JCheckBox box : catalogCheckboxes) {
                    box.setSelected(true);
                }
            }
        });

        // Deselect All listener
        checkboxdeall.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (JCheckBox box : catalogCheckboxes) {
                    box.setSelected(false);
                }
            }
        });

        // Creates a "WiseView" checkbox area
        JCheckBox checkboxwise = new JCheckBox("WiseView");
        checkboxwise.setFont(new Font("Times New Roman", Font.PLAIN, 9));
        checkboxwise.setBounds(2, 175, 80, 20);
        mainControlPanel.add(checkboxwise);

        // Creates a field for the button "Run"
        JButton runbutton = new JButton("Run");
        runbutton.setBounds(75, 165, 200, 30);
        mainControlPanel.add(runbutton);
        // ------------------------------------------------------------- //


        
        // Listens for when the "Run" button is clicked
        // ------------------------------------------------------------- //
        runbutton.addActionListener(new ActionListener() {
            private static int fileCount = 1;
            @Override
            public void actionPerformed(ActionEvent e) {
                runbutton.setEnabled(false);

                // Gets the RA, DEC, and RADIUS inputs
                String raText = raField.getText();
                String decText = decField.getText();
                String radiusText = radiusField.getText();

                // Verifies that the RA, DEC, and RADIUS inputs are flaots
                try {
                    // Checks if "WiseView" checkbox is on
                    boolean wisevar = checkboxwise.isSelected();

                    // Verify the RADIUS are within the given bounds
                    Float radiusFloat = Float.parseFloat(radiusText);
    
                    // Adds all checked catalogs to a master list
                    List<String> selectedCheckboxes = new ArrayList<>();
                    if (checkBox1.isSelected()) selectedCheckboxes.add("CatWISE");
                    if (checkBox2.isSelected()) selectedCheckboxes.add("AllWISE");
                    if (checkBox3.isSelected()) selectedCheckboxes.add("Gaia");
                    if (checkBox4.isSelected()) selectedCheckboxes.add("VHS");
                    if (checkBox5.isSelected()) selectedCheckboxes.add("UKIDSS");
                    if (checkBox6.isSelected()) selectedCheckboxes.add("2MASS");
                    if (checkBox7.isSelected()) selectedCheckboxes.add("SDSS");
                    if (checkBox8.isSelected()) selectedCheckboxes.add("PanSTARRS");
                    if (checkBox9.isSelected()) selectedCheckboxes.add("NSC");
                    if (checkBox10.isSelected()) selectedCheckboxes.add("GALEX");

                    // Checks if any catalog is selected
                    if (selectedCheckboxes.isEmpty()) {
                        JOptionPane.showMessageDialog(runbutton, "Select a Catalog!", "Input Error", JOptionPane.ERROR_MESSAGE);
                    } else {

                        boolean multivar = checkboxmultiquery.isSelected();

                        // Use a relative path from the location of the JAR
                        File scriptFile = new File(appDir, "/resources/WRAP.py");
                        File pythonExec = new File(appDir, "/resources/myenv/bin/python3.8");

                        // Add the command and needed information for "WRAP.py"
                        List<String> command = new ArrayList<>();
                        command.add(pythonExec.getAbsolutePath());
                        command.add(scriptFile.getAbsolutePath());
                        String csvpath = "";
                        if (multivar){
                            csvpath = fileField.getText();
                            command.add(String.valueOf(csvpath));
                            command.add(String.valueOf(0));
                        }else{
                            Float raFloat = Float.parseFloat(raText);
                            Float decFloat = Float.parseFloat(decText);
                            if (raFloat <= 360 && raFloat >= 0 && decFloat >= -90 && decFloat <= 90 && radiusFloat >= 50 && radiusFloat <= 500){
                                command.add(String.valueOf(raFloat));
                                command.add(String.valueOf(decFloat));
                            } else{
                                JOptionPane.showMessageDialog(runbutton, "Input Correct RA (BETWEEN 0 AND 360)/Decl (BETWEEN -90 AND 90)/Radius (BETWEEN 50 AND 500)", "Input Error", JOptionPane.ERROR_MESSAGE);
                            }
                        }
                        command.add(String.valueOf(radiusFloat));
                        command.add(String.valueOf(multivar));
                        command.add(String.valueOf(wisevar));
                        command.addAll(selectedCheckboxes);

                        

                        // Start the process
                        ProcessBuilder pb = new ProcessBuilder(command);
                        pb.redirectErrorStream(true);
                        Process process = pb.start();

                        // Read output
                        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        StringBuilder output = new StringBuilder();
                        String line;

                        while ((line = reader.readLine()) != null) {
                            output.append(line);
                        }

                        // Store the result in a variable
                        String resultJson = output.toString();
                        String filePrefix = "WRAP_Query_";
                        String fileName = filePrefix + fileCount;
                        
                        // Stores user click data into left panel for tables
                        fileLoaderPanel.loadJsonData(resultJson, fileName, multivar, csvpath);
                        fileCount++;
                    }

                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(runbutton, ex.getMessage(), "Input Error", JOptionPane.ERROR_MESSAGE);
                }
                runbutton.setEnabled(true);
            }
        });
        // ------------------------------------------------------------- //



        // Creates a "?" button that pops up with markdown code for help
        // ------------------------------------------------------------- //
        JButton help_button = new JButton("?");
        help_button.setBounds(315, 170, 20, 20);
        mainControlPanel.add(help_button);
        help_button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    // Path to your .md file
                    Path readmePath = Paths.get(appDir + "/resources/wrap_help.md");
                    String content = Files.readString(readmePath);

                    // Remove unwanted symbols or characters
                    content = content.replace("java", ""); 
                    content = content.replaceAll("[^\\x00-\\x7F]", "");

                    // Display the cleaned content
                    JTextArea textArea = new JTextArea(content);
                    textArea.setEditable(false);
                    textArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
                    JScrollPane scrollPane = new JScrollPane(textArea);
                    scrollPane.setPreferredSize(new Dimension(500, 400));

                    JOptionPane.showMessageDialog(frame, scrollPane, "WRAP Help", JOptionPane.INFORMATION_MESSAGE);
                } catch (IOException ex) {
                    JOptionPane.showMessageDialog(frame, "Error loading help file: " + ex.getMessage(),
                            "Help Load Error", JOptionPane.ERROR_MESSAGE);
                }
            }
        });
        // ------------------------------------------------------------- //

        // Adds all right side information to the window and shows the window
        frame.add(mainControlPanel);
        frame.setVisible(true);
    }
// ------------------------------------------------------------- //


    // FileLoaderPanel class for file loading and displaying file names
    // ------------------------------------------------------------- //
    static class FileLoaderPanel extends JPanel {
        // Sets up the needed information
        // ------------------------------------------------------------- //
        private JList<String> fileList;
        private DefaultListModel<String> listModel;
        private JButton loadButton;
        private JButton removeButton;
        private JButton saveButton;

        private Map<String, Map<String, String>> tableData = new HashMap<>();
        // ------------------------------------------------------------- //



        // Creates panel to store all loaded data
        // ------------------------------------------------------------- //
        public FileLoaderPanel() {
            setLayout(new BorderLayout());

            // Creates the window and allows for scrolling
            listModel = new DefaultListModel<>();
            fileList = new JList<>(listModel);
            fileList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
            JScrollPane scrollPane = new JScrollPane(fileList);
            scrollPane.setBorder(BorderFactory.createTitledBorder("Table List"));

            // Adds the "Save" and "Remove" buttons
            JPanel buttonPanel = new JPanel();
            saveButton = new JButton("Save");
            removeButton = new JButton("Remove");

            // Makes these button look pretty
            removeButton.setFont(new Font("Times New Roman", Font.PLAIN, 12));
            saveButton.setFont(new Font("Times New Roman", Font.PLAIN, 12));
            removeButton.setPreferredSize(new Dimension(64, 30));
            saveButton.setPreferredSize(new Dimension(64, 30));

            // Add the above buttons to the panel
            buttonPanel.add(removeButton);
            buttonPanel.add(saveButton);
            add(scrollPane, BorderLayout.CENTER);
            add(buttonPanel, BorderLayout.SOUTH);

            // Listens for when either button is pressed
            removeButton.addActionListener(e -> removeFile());
            saveButton.addActionListener(e -> saveFile());
        }
        // ------------------------------------------------------------- //



        // Created to show the loaded data files from the user click data
        // ------------------------------------------------------------- //
        public void loadJsonData(String jsonString, String tableName, Boolean multivar, String csvpath) {
            Map<String, String> tableMap = new LinkedHashMap<>();
            jsonString = jsonString.trim();

            // Find the start of the actual data
            int startIndex = jsonString.indexOf("input_ra");
            // if (startIndex == -1) {
            //     JOptionPane.showMessageDialog(this, "No Data Available \n (e.g. Servers Are Currently Down)", "Load Error", JOptionPane.ERROR_MESSAGE);
            //     return;
            // }
        
            // Trim off any log junk before 'input_ra'
            jsonString = jsonString.substring(startIndex);
        
            // Split the clean string
            String[] parts = jsonString.split(",");

            int rowCount = 0;
            if (multivar){
                try (BufferedReader br = new BufferedReader(new FileReader(csvpath))) {
                    while (br.readLine() != null) {
                        rowCount++;
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }else{
                rowCount = 2;
            }
        
            int total = parts.length;
            if (total % rowCount != 0) {
                JOptionPane.showMessageDialog(this, "Data format error: uneven number of columns and values", "Load Error", JOptionPane.ERROR_MESSAGE);
                return;
            }
        
            int columns = total / rowCount;
            String[] keys = new String[columns];

            // Extract keys from the first row
            for (int i = 0; i < columns; i++) {
                keys[i] = parts[i].trim();
            }

            // Loop through the remaining rows and fill in values
            for (int row = 1; row < rowCount; row++) {
                for (int col = 0; col < columns; col++) {
                    String key = keys[col] + "_row" + (row - 1);
                    String value = parts[row * columns + col].trim();
                    tableMap.put(key, value);
                }
            }

            tableData.put(tableName, tableMap);
            listModel.addElement(tableName);
        }      
        // ------------------------------------------------------------- //         



        // Removes unwanted files when user clicked
        // ------------------------------------------------------------- //
        private void removeFile() {
            int selectedIndex = fileList.getSelectedIndex();
            if (selectedIndex != -1) {
                String key = fileList.getSelectedValue();
                listModel.remove(selectedIndex);
                tableData.remove(key);
            } else {
                JOptionPane.showMessageDialog(this, "No File Selected", "Removal Error", JOptionPane.ERROR_MESSAGE);
            }
        }
        // ------------------------------------------------------------- //



        // When user clicks "Save" it provides a popup to save user click data to a csv
        // ------------------------------------------------------------- //
        private void saveFile() {
            int selectedIndex = fileList.getSelectedIndex();
            if (selectedIndex != -1) {
                String tableName = fileList.getSelectedValue();
                Map<String, String> table = tableData.get(tableName);

                // Reorganize keys by base name and row index
                Map<Integer, Map<String, String>> rows = new TreeMap<>();
                Set<String> baseKeys = new LinkedHashSet<>();

                for (String fullKey : table.keySet()) {
                    String[] parts = fullKey.split("_row");
                    String baseKey = parts[0];
                    int rowIndex = (parts.length > 1) ? Integer.parseInt(parts[1]) : 0;

                    baseKeys.add(baseKey);
                    rows.computeIfAbsent(rowIndex, k -> new LinkedHashMap<>()).put(baseKey, table.get(fullKey));
                }

                JFileChooser fileChooser = new JFileChooser();
                fileChooser.setSelectedFile(new File(tableName + ".csv"));
                int returnValue = fileChooser.showSaveDialog(this);

                if (returnValue == JFileChooser.APPROVE_OPTION) {
                    File saveFile = fileChooser.getSelectedFile();
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(saveFile))) {
                        // Write header
                        writer.write(String.join(",", baseKeys));
                        writer.newLine();

                        // Write rows
                        for (Map<String, String> row : rows.values()) {
                            List<String> values = new ArrayList<>();
                            for (String key : baseKeys) {
                                values.add(row.getOrDefault(key, ""));
                            }
                            writer.write(String.join(",", values));
                            writer.newLine();
                        }

                        JOptionPane.showMessageDialog(this, "File saved successfully.", "Success", JOptionPane.INFORMATION_MESSAGE);
                    } catch (IOException ex) {
                        JOptionPane.showMessageDialog(this, "Error saving file: " + ex.getMessage(), "File Save Error", JOptionPane.ERROR_MESSAGE);
                    }
                }
            } else {
                JOptionPane.showMessageDialog(this, "No File Selected", "Save Error", JOptionPane.ERROR_MESSAGE);
            }
        }
        // ------------------------------------------------------------- //
    }
}
