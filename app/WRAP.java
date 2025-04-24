import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
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
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
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
import javax.swing.SwingConstants;


public class WRAP {
    public static void main(String[] args) {
        JFrame frame = new JFrame("WRAP");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(590, 255);
        frame.setLayout(null);
        frame.setResizable(false);

        String jarPath = "";
        try {
            jarPath = new File(WRAP.class.getProtectionDomain().getCodeSource().getLocation().toURI()).getPath();
        } catch (URISyntaxException e) {
            e.printStackTrace(); // or handle more gracefully in a GUI
        }
        File appDir = new File(jarPath).getParentFile(); // Points to "app/"

        // File Loader Panel
        FileLoaderPanel fileLoaderPanel = new FileLoaderPanel();
        fileLoaderPanel.setBounds(10, 10, 220, 220); 
        frame.add(fileLoaderPanel);

        // WRAP Data Panel
        JPanel mainControlPanel = new JPanel();
        mainControlPanel.setLayout(null);
        mainControlPanel.setPreferredSize(new Dimension(200, 350));
        mainControlPanel.setBounds(235, 10, 345, 200);
        mainControlPanel.setBorder(BorderFactory.createTitledBorder("Input Requirements"));

        // RA, DEC, RADIUS input
        JPanel coordsPanel = new JPanel(new GridLayout(1, 4));

        JTextField raField = new JTextField("R.A.", 20);
        raField.setForeground(Color.GRAY);
        raField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (raField.getText().equals("R.A.")) {
                    raField.setText("");
                    raField.setForeground(Color.BLACK);
                }
            }

            @Override
            public void focusLost(FocusEvent e) {
                if (raField.getText().isEmpty()) {
                    raField.setText("R.A.");
                    raField.setForeground(Color.GRAY);
                }
            }
        });


        JTextField decField = new JTextField("Decl.", 20);
        decField.setForeground(Color.GRAY);
        decField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (decField.getText().equals("Decl.")) {
                    decField.setText("");
                    decField.setForeground(Color.BLACK);
                }
            }

            @Override
            public void focusLost(FocusEvent e) {
                if (decField.getText().isEmpty()) {
                    decField.setText("Decl.");
                    decField.setForeground(Color.GRAY);
                }
            }
        });

        JTextField radiusField = new JTextField("Search Radius", 20);
        radiusField.setForeground(Color.GRAY);
        radiusField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (radiusField.getText().equals("Search Radius")) {
                    radiusField.setText("");
                    radiusField.setForeground(Color.BLACK);
                }
            }

            @Override
            public void focusLost(FocusEvent e) {
                if (radiusField.getText().isEmpty()) {
                    radiusField.setText("Search Radius");
                    radiusField.setForeground(Color.GRAY);
                }
            }
        });
        
        coordsPanel.add(raField);
        coordsPanel.add(decField);
        coordsPanel.add(radiusField);
        coordsPanel.setBounds(5, 20, 335, 25);
        mainControlPanel.add(coordsPanel);

        // Catalog Panel with checkboxes
        JPanel catalogPanel = new JPanel(new GridLayout(3, 3));
        JCheckBox checkBox1 = new JCheckBox("CatWISE");
        JCheckBox checkBox2 = new JCheckBox("AllWISE");
        JCheckBox checkBox3 = new JCheckBox("Gaia");
        JCheckBox checkBox4 = new JCheckBox("VISTA");
        JCheckBox checkBox5 = new JCheckBox("WFCAM");
        JCheckBox checkBox6 = new JCheckBox("2MASS");
        JCheckBox checkBox7 = new JCheckBox("PanSTARRS");
        JCheckBox checkBox8 = new JCheckBox("NSC");
        JCheckBox checkBox9 = new JCheckBox("GALEX");

        // Set horizontal alignment for the text in the checkboxes
        checkBox1.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox1.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox2.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox2.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox3.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox3.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox4.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox4.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox5.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox5.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox6.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox6.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox7.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox7.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox8.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox8.setFont(new Font("Times New Roman", Font.PLAIN, 15));
        checkBox9.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox9.setFont(new Font("Times New Roman", Font.PLAIN, 15));

        // Add checkboxes to the catalog panel
        catalogPanel.add(checkBox1);
        catalogPanel.add(checkBox2);
        catalogPanel.add(checkBox3);
        catalogPanel.add(checkBox4);
        catalogPanel.add(checkBox5);
        catalogPanel.add(checkBox6);
        catalogPanel.add(checkBox7);
        catalogPanel.add(checkBox8);
        catalogPanel.add(checkBox9);

        // Set the bounds of the catalog panel and add it to the frame
        catalogPanel.setBounds(3, 40, 335, 100);
        mainControlPanel.add(catalogPanel);

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
        optionsPanel.setBounds(100, 135, 150, 20);
        mainControlPanel.add(optionsPanel);

        // Array or list to hold all catalog checkboxes
        JCheckBox[] catalogCheckboxes = {
            checkBox1, checkBox2, checkBox3, checkBox4, checkBox5,
            checkBox6, checkBox7, checkBox8, checkBox9
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

        JCheckBox checkboxwise = new JCheckBox("WiseView");
        checkboxwise.setFont(new Font("Times New Roman", Font.PLAIN, 10));
        checkboxwise.setBounds(2, 170, 80, 20);
        mainControlPanel.add(checkboxwise);

        // Buttons
        JButton runbutton = new JButton("Run");
        runbutton.setBounds(70, 160, 200, 30);
        mainControlPanel.add(runbutton);
        
        runbutton.addActionListener(new ActionListener() {
            private static int fileCount = 1;
            @Override
            public void actionPerformed(ActionEvent e) {
                runbutton.setEnabled(false);
                String raText = raField.getText();
                String decText = decField.getText();
                String radiusText = radiusField.getText();
                try {
                    Float raFloat = Float.parseFloat(raText);
                    Float decFloat = Float.parseFloat(decText);
                    Float radiusFloat = Float.parseFloat(radiusText);

                    if (checkboxwise.isSelected()){
                        try {
                            URI uri = new URI("http://byw.tools/wiseview#ra=" + raText + "&dec=" + decText + "&size=" + radiusText + "&band=3&speed=164.06&minbright=-10.0000&maxbright=80&window=0.75&diff_window=1&linear=1&color=&zoom=8.5&border=1&gaia=1&invert=1&maxdyr=0&scandir=0&neowise=0&diff=0&outer_epochs=0&unique_window=1&smooth_scan=0&shift=0&pmra=0&pmdec=0&synth_a=0&synth_a_sub=0&synth_a_ra=&synth_a_dec=&synth_a_w1=&synth_a_w2=&synth_a_pmra=0&synth_a_pmdec=0&synth_a_mjd=&synth_b=0&synth_b_sub=0&synth_b_ra=&synth_b_dec=&synth_b_w1=&synth_b_w2=&synth_b_pmra=0&synth_b_pmdec=0&synth_b_mjd=&smdet_coadd_id=1863p620&smdet_mask_idx=3&smdet_obj_px=&smdet_obj_py=");

                            if (Desktop.isDesktopSupported()) {
                                Desktop desktop = Desktop.getDesktop();
                                desktop.browse(uri);
                            }
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }

                    List<String> selectedCheckboxes = new ArrayList<>();
                    if (checkBox1.isSelected()) selectedCheckboxes.add("CatWISE");
                    if (checkBox2.isSelected()) selectedCheckboxes.add("AllWISE");
                    if (checkBox3.isSelected()) selectedCheckboxes.add("Gaia");
                    if (checkBox4.isSelected()) selectedCheckboxes.add("VISTA");
                    if (checkBox5.isSelected()) selectedCheckboxes.add("WFCAM");
                    if (checkBox6.isSelected()) selectedCheckboxes.add("2MASS");
                    if (checkBox7.isSelected()) selectedCheckboxes.add("PanSTARRS");
                    if (checkBox8.isSelected()) selectedCheckboxes.add("NSC");
                    if (checkBox9.isSelected()) selectedCheckboxes.add("GALEX");

                    // Use a relative path from the location of the JAR
                    File scriptFile = new File(appDir, "/resources/WRAP.py");
                    File pythonExec = new File(appDir, "/resources/wrapenv/bin/python3.8");

                    List<String> command = new ArrayList<>();
                    command.add(pythonExec.getAbsolutePath());
                    command.add(scriptFile.getAbsolutePath());
                    command.add(String.valueOf(raFloat));
                    command.add(String.valueOf(decFloat));
                    command.add(String.valueOf(radiusFloat));
                    command.addAll(selectedCheckboxes);

                    // Start the process
                    ProcessBuilder pb = new ProcessBuilder(command);
                    pb.redirectErrorStream(true); // Merge stderr into stdout
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
                    fileLoaderPanel.loadJsonData(resultJson, fileName);
                    fileCount++;
                    
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(runbutton, "Error: " + ex.getMessage(), "Input Error", JOptionPane.ERROR_MESSAGE);
                }
                runbutton.setEnabled(true);
            }
        });


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

                    // Optional: sanitize content (e.g., remove unwanted symbols or characters)
                    content = content.replace("java", ""); // Removes the "java" symbol if it appears.
                    content = content.replaceAll("[^\\x00-\\x7F]", ""); // Removes non-ASCII characters (if necessary)

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

        frame.add(mainControlPanel);
        frame.setVisible(true);
    }

    // FileLoaderPanel class for file loading and displaying file names
    static class FileLoaderPanel extends JPanel {
        private JList<String> fileList;
        private DefaultListModel<String> listModel;
        private JButton loadButton;
        private JButton removeButton;
        private JButton saveButton;

        private Map<String, Map<String, String>> tableData = new HashMap<>();

        public FileLoaderPanel() {
            setLayout(new BorderLayout());

            listModel = new DefaultListModel<>();
            fileList = new JList<>(listModel);
            fileList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
            JScrollPane scrollPane = new JScrollPane(fileList);
            scrollPane.setBorder(BorderFactory.createTitledBorder("Table List"));

            JPanel buttonPanel = new JPanel();
            saveButton = new JButton("Save");
            removeButton = new JButton("Remove");

            removeButton.setFont(new Font("Times New Roman", Font.PLAIN, 12));
            saveButton.setFont(new Font("Times New Roman", Font.PLAIN, 12));

            removeButton.setPreferredSize(new Dimension(64, 30));
            saveButton.setPreferredSize(new Dimension(64, 30));

            buttonPanel.add(removeButton);
            buttonPanel.add(saveButton);

            add(scrollPane, BorderLayout.CENTER);
            add(buttonPanel, BorderLayout.SOUTH);

            removeButton.addActionListener(e -> removeFile());
            saveButton.addActionListener(e -> saveFile());
        }

        public void loadJsonData(String jsonString, String tableName) {
            Map<String, String> tableMap = new LinkedHashMap<>();
            jsonString = jsonString.trim();

            if (jsonString.startsWith("[") && jsonString.endsWith("]")) {
                jsonString = jsonString.substring(1, jsonString.length() - 1);
            }
            if (jsonString.startsWith("{") && jsonString.endsWith("}")) {
                jsonString = jsonString.substring(1, jsonString.length() - 1);
            }

            String[] entries = jsonString.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)");
            for (String entry : entries) {
                String[] pair = entry.split(":", 2);
                if (pair.length == 2) {
                    String key = pair[0].trim().replaceAll("^\"|\"$", "");
                    String value = pair[1].trim().replaceAll("^\"|\"$", "");
                    tableMap.put(key, value);
                }
            }

            tableData.put(tableName, tableMap);
            listModel.addElement(tableName);
        }

        private void loadFiles() {
            JFileChooser fileChooser = new JFileChooser();
            fileChooser.setMultiSelectionEnabled(true);
            int returnValue = fileChooser.showOpenDialog(this);

            if (returnValue == JFileChooser.APPROVE_OPTION) {
                File[] files = fileChooser.getSelectedFiles();
                for (File file : files) {
                    listModel.addElement(file.getName());
                }
            }
        }

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

        private void saveFile() {
            int selectedIndex = fileList.getSelectedIndex();
            if (selectedIndex != -1) {
                String tableName = fileList.getSelectedValue();
                Map<String, String> table = tableData.get(tableName);

                JFileChooser fileChooser = new JFileChooser();
                fileChooser.setSelectedFile(new File(tableName + ".csv"));
                int returnValue = fileChooser.showSaveDialog(this);

                if (returnValue == JFileChooser.APPROVE_OPTION) {
                    File saveFile = fileChooser.getSelectedFile();
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(saveFile))) {
                        for (String key : table.keySet()) {
                            writer.write(key + ",");
                        }
                        writer.newLine();

                        for (String value : table.values()) {
                            writer.write(value + ",");
                        }
                        writer.newLine();

                        JOptionPane.showMessageDialog(this, "File saved successfully.", "Success", JOptionPane.INFORMATION_MESSAGE);
                    } catch (IOException ex) {
                        JOptionPane.showMessageDialog(this, "Error saving file: " + ex.getMessage(), "File Save Error", JOptionPane.ERROR_MESSAGE);
                    }
                }
            } else {
                JOptionPane.showMessageDialog(this, "No File Selected", "Save Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }
}
