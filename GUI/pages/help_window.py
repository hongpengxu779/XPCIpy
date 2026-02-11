import tkinter as tk
from tkinter import ttk
import GUI.i18n as i18n

class HelpWindow:
    
    def __init__(self, master):
        self.master = master
        self.build_window()

    def build_window(self):
        self.win = tk.Toplevel(self.master)
        self.win.title(i18n.HELP_TITLE)
        self.win.geometry("800x600")

        notebook = ttk.Notebook(self.win)
        notebook.pack(fill="both", expand=True)

        self.tab_sim = ttk.Frame(notebook)
        self.tab_rec = ttk.Frame(notebook)
        self.tab_rec_batch = ttk.Frame(notebook)
        self.tab_about = ttk.Frame(notebook)
        self.tab_license = ttk.Frame(notebook)
        

        notebook.add(self.tab_sim, text=i18n.HELP_TAB_SIM)
        notebook.add(self.tab_rec, text=i18n.HELP_TAB_REC)
        notebook.add(self.tab_rec_batch, text=i18n.HELP_TAB_BATCH)
        notebook.add(self.tab_about, text=i18n.HELP_TAB_ABOUT)

        # Fill content
        self._populate_simulation_guide()
        self._populate_reconstruction_guide()
        self._populate_reconstruction_batch()
        self._populate_about_tab()
  
    def _make_text_widget(self, parent):
        text = tk.Text(parent, wrap="word", bg="#222222", fg="white")
        text.pack(fill="both", expand=True, padx=10, pady=10)

        text.tag_configure("title", font=("Arial", 14, "bold"), foreground="#ffd966")
        text.tag_configure("subtitle", font=("Arial", 12, "bold"), foreground="#ffffff")
        text.tag_configure("bullet", lmargin1=25, lmargin2=40)

        return text

    def _populate_simulation_guide(self):
        text = self._make_text_widget(self.tab_sim)

        text.insert("end", i18n.HELP_SIM_TITLE + "\n", "title")

        text.insert("end", "\n" + i18n.HELP_INLINE_TITLE + "\n", "subtitle")
        text.insert("end", i18n.HELP_INLINE_DESC)
        
        for b in i18n.HELP_INLINE_BULLETS:
            text.insert("end", b, "bullet")

        text.insert("end", "\n" + i18n.HELP_TL_SIM_TITLE + "\n", "subtitle")
        text.insert("end", i18n.HELP_TL_SIM_DESC)
        
        for b in i18n.HELP_TL_SIM_BULLETS:
            text.insert("end", b, "bullet")

        text.insert("end", "\n" + i18n.HELP_CHECK_TL_TITLE + "\n", "subtitle")
        text.insert("end", i18n.HELP_CHECK_TL_DESC)

        text.config(state="disabled")

    def _populate_reconstruction_guide(self):
        text = self._make_text_widget(self.tab_rec)

        text.insert("end", i18n.HELP_REC_TITLE + "\n", "title")

        text.insert("end", "\n" + i18n.HELP_REC_WORKFLOW + "\n", "subtitle")
        text.insert("end", i18n.HELP_REC_DESC)
        for b in i18n.HELP_REC_BULLETS:
            text.insert("end", b, "bullet")

        text.insert("end", "\n" + i18n.HELP_REC_TOOLS_TITLE + "\n", "subtitle")
        text.insert("end", i18n.HELP_REC_TOOLS_DESC)
        for b in i18n.HELP_REC_TOOLS_BULLETS:
            text.insert("end", b, "bullet")

        text.config(state="disabled")
        
    def _populate_reconstruction_batch(self):
        text = self._make_text_widget(self.tab_rec_batch)
        text.insert("end", "\n" + i18n.HELP_BATCH_TITLE + "\n", "title")

        text.insert("end", "\n" + i18n.HELP_BATCH_OVERVIEW + "\n", "subtitle")
        text.insert("end", i18n.HELP_BATCH_OVERVIEW_DESC)

        text.insert("end", "\n" + i18n.HELP_BATCH_FOLDER_TITLE + "\n", "subtitle")
        text.insert("end", i18n.HELP_BATCH_FOLDER_DESC)
        for b in i18n.HELP_BATCH_FOLDER_BULLETS:
            text.insert("end", b, "bullet")
        text.insert("end", i18n.HELP_BATCH_FOLDER_NOTE)

        text.insert("end", "\n" + i18n.HELP_BATCH_WHAT_TITLE + "\n", "subtitle")
        text.insert("end", i18n.HELP_BATCH_WHAT_DESC)
        for b in i18n.HELP_BATCH_WHAT_BULLETS:
            text.insert("end", b, "bullet")

        text.insert("end", "\n" + i18n.HELP_BATCH_OUTPUT_TITLE + "\n", "subtitle")
        for b in i18n.HELP_BATCH_OUTPUT_BULLETS:
            text.insert("end", b, "bullet")

        text.insert("end", "\n" + i18n.HELP_BATCH_HOW_TITLE + "\n", "subtitle")
        for b in i18n.HELP_BATCH_HOW_BULLETS:
            text.insert("end", b, "bullet")
        text.insert("end", i18n.HELP_BATCH_HOW_NOTE)

        text.insert("end", "\n" + i18n.HELP_BATCH_NOTES_TITLE + "\n", "subtitle")
        for b in i18n.HELP_BATCH_NOTES_BULLETS:
            text.insert("end", b, "bullet")
        
        text.config(state="disabled")

    def _populate_about_tab(self):
        text = self._make_text_widget(self.tab_about)

        text.insert("end", i18n.HELP_ABOUT_TITLE + "\n", "title")
        text.insert("end", i18n.HELP_ABOUT_DESC)

        text.insert("end", "\n" + i18n.HELP_ABOUT_MODULES + "\n", "subtitle")
        text.insert("end", i18n.HELP_ABOUT_MODULES_DESC)

        text.insert("end", "\n" + i18n.HELP_ABOUT_AUTHOR + "\n", "subtitle")
        text.insert("end", "  Víctor Sánchez Lara\n")
        text.insert("end", "  Email: vicsan05@ucm.es\n")
    
        text.insert("end", "\n" + i18n.HELP_ABOUT_PROJECT + "\n", "subtitle")
        text.insert("end", i18n.HELP_ABOUT_PROJECT_DESC)

        text.config(state="disabled")

