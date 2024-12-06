import pandas as pd
import argparse
import argcomplete
import csv
import os
import fpdf
from fpdf import FPDF
import pdf2image
from pdf2image import convert_from_path
from datetime import datetime
import reportlab
import pdfrw
from pdfrw import PdfReader, PdfWriter, PdfName
#for fillable boxes
from reportlab.pdfgen import canvas  
from reportlab.pdfbase import pdfform  
from reportlab.lib.colors import magenta, pink, blue, green
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfbase.pdfform import textFieldAbsolute
from reportlab.platypus import PageBreak
from reportlab.platypus import SimpleDocTemplate

import pdf_highlighter as pdf_highlighter

#descriptions and argument assignment
parser = argparse.ArgumentParser(description='details',
         usage='use "%(prog)s --help" for more information',
         formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--info', default=None,
                    help='''
                    For use with redmine automator to create GROBI reports
                    ''')
parser.add_argument('-s', dest='seqid', type=str, required=True, help="""SEQ-ID""")
parser.add_argument('-i_rmlst', dest='rmlst_infile', type=str, required=True, help="""rmlst_file.txt""")
parser.add_argument('-i_anib', dest='ANIb_infile', type=str, required=True, help="""ANIb_input.tab""")
parser.add_argument('-anib_mat', dest='ANIb_matrix', type=str, required=True, help="""ANIb_percentage_identity.png""")
parser.add_argument('-i_anim', dest='ANIm_infile', type=str, required=True, help="""ANIm_input.tab""")
parser.add_argument('-anim_mat', dest='ANIm_matrix', type=str, required=True, help="""ANIm_percentage_identity.png""")
parser.add_argument('-i_mash', dest='MASH_infile', type=str, required=True, help="""MASH-distances.tab""")
parser.add_argument('-drep_tree', dest='dRep_tree', type=str, required=True, help="""Primary_clustering_dendrogram.pdf""")
parser.add_argument('-o', dest='outfile', type=str, required=True, help="""output.pdf""")
argcomplete.autocomplete(parser)
args = parser.parse_args()

#parse the pubmlst outfile
rank_line = [0]
taxon_line = [1]
support_line = [2]
taxonomy_line = [3]
o_file = open(args.rmlst_infile)
for index, line in enumerate(o_file):
    if index in rank_line:
        Rank = line.split(':')[1]
    if index in taxon_line:
        Taxon = line.split(':')[1]
    if index in support_line:
        Support = line.split(':')[1]
        Supp_num = Support.rstrip("%") #TODO: come back to this and remove percent sign and convert to a number?
    if index in taxonomy_line:
        Taxonomy = line.split('axonomy')[1]

#parse the pyani outfiles
anib_df = pd.read_csv(args.ANIb_infile, sep='\t')
ANIb_referencerow = anib_df.sort_values(by=['{}'.format(args.seqid)], ascending=False).head(2) #we want the second row, as the first should be the reference against itself (unless something has gone horribly wrong)
ANIb_reference_id = ANIb_referencerow.iloc[1]['Unnamed: 0']
ANIb_ANIvalue = ANIb_referencerow.iloc[1]['{}'.format(args.seqid)]

anim_df = pd.read_csv(args.ANIm_infile, sep='\t')
ANIm_referencerow = anim_df.sort_values(by=['{}'.format(args.seqid)], ascending=False).head(2) #we want the second row, as the first should be the reference against itself
ANIm_reference_id = ANIm_referencerow.iloc[1]['Unnamed: 0']
ANIm_ANIvalue = ANIm_referencerow.iloc[1]['{}'.format(args.seqid)]

#parse the mash outfile
mash_df = pd.read_csv("{}".format(args.MASH_infile), sep='\t', names=["Reference-ID","Query-ID","MASH-Distance","P-value","Matching-hashes"])
referencerow = mash_df.sort_values(by=['MASH-Distance'], ascending=True).head(1)
referenceid = referencerow.iloc[0]['Reference-ID']
refidtable = os.path.basename(referenceid).split(".")[0]
minhashes = referencerow.iloc[0]['Matching-hashes']

#highlight the sequence on the MASH tree using the pdf_highlighter script
pdf_highlighter.process_data('{}'.format(args.dRep_tree),'Primary_clustering_dendrogram_highlighted.pdf',"{}".format(args.seqid))
#also add a red box around the sequenceID on the MASH tree, because why not?!
pdf_highlighter.process_data('Primary_clustering_dendrogram_highlighted.pdf','Primary_clustering_dendrogram_highlighted.pdf',"{}".format(args.seqid),None,'Frame')

#convert the modified dRep primary clustering tree to an image file
images = convert_from_path('Primary_clustering_dendrogram_highlighted.pdf')
for i in range(len(images)):
    #save the pages as images in the pdf
    images[i].save('Primary_clustering_dendrogram'+ str(i) +'.png','PNG')

#create the output summary table as a pdf
# Margin
m = 10
# Page width: Width of A4 is 210mm
pw = 210 - 2*m
# cell height
ch = 8

pdf = fpdf.FPDF(orientation="portrait", format=(210,80))
pdf.add_page()
#table with each element individually set.. maybe this way can add "if" statements for values?
pdf.set_font('Arial', 'B', 16)
pdf.cell(w=30, h=ch, txt="Analysis", border=1, ln=0, align='C')
pdf.cell(w=130, h=ch, txt="Matching Reference", border=1, ln=0, align='C')
pdf.cell(w=32, h=ch, txt="Support", border=1, ln=1, align='C')

#now put the values in
pdf.set_font('Arial', '', 9)
#rMLST
pdf.cell(w=30, h=ch, txt="rMLST", border=1, ln=0, align='C')
pdf.cell(w=130, h=ch, txt="{}".format(Taxon), border=1, ln=0, align='C')
pdf.cell(w=32, h=ch, txt="{}".format(Support), border=1, ln=1, align='C')

#ANIb
pdf.cell(w=30, h=ch, txt="ANIb", border=1, ln=0, align='C')
if ANIb_ANIvalue < 0.95:
    pdf.set_fill_color(r=203,g=91,b=76)
    pdf.cell(w=130, h=ch, txt="{}".format(ANIb_reference_id), border=1, ln=0, align='C')
    pdf.cell(w=32, h=ch, txt="{}".format(ANIb_ANIvalue), border=1, ln=1, align='C', fill=True)
else:
    pdf.cell(w=130, h=ch, txt="{}".format(ANIb_reference_id), border=1, ln=0, align='C')
    pdf.cell(w=32, h=ch, txt="{}".format(ANIb_ANIvalue), border=1, ln=1, align='C')

#ANIm
pdf.cell(w=30, h=ch, txt="ANIm", border=1, ln=0, align='C')
if ANIm_ANIvalue < 0.95:
    pdf.set_fill_color(r=203,g=91,b=76)
    pdf.cell(w=130, h=ch, txt="{}".format(ANIm_reference_id), border=1, ln=0, align='C')
    pdf.cell(w=32, h=ch, txt="{}".format(ANIm_ANIvalue), border=1, ln=1, align='C', fill=True)
else:
    pdf.cell(w=130, h=ch, txt="{}".format(ANIb_reference_id), border=1, ln=0, align='C')
    pdf.cell(w=32, h=ch, txt="{}".format(ANIm_ANIvalue), border=1, ln=1, align='C')

#MASH
pdf.cell(w=30, h=ch, txt="MASH", border=1, ln=0, align='C')
pdf.cell(w=130, h=ch, txt="{}".format(refidtable), border=1, ln=0, align='C')
pdf.cell(w=32, h=ch, txt="{}".format(minhashes), border=1, ln=1, align='C')

pdf.output('GROBI_table.pdf', 'F')

#convert pdf to image
images = convert_from_path('GROBI_table.pdf')
for i in range(len(images)):
    #save the pages as images in the pdf
    #images[i].save('page'+ str(i) +'.jpg','JPEG')
    images[i].save('GROBI_table'+ '.png','PNG')

#now create the ANI figure pages
# Margin
m = 10
# Page width: Width of A4 is 210mm
pw = 210 - 2*m
# cell height
ch = 8
class PDF(FPDF):
    def __init__(self):
        super().__init__()
    #def header(self):
    #    self.set_font('Arial', '', 12)
    #    self.cell(0, 8, 'Header', 0, 1, 'C')
#    def footer(self):
#        self.set_y(-15)
#        self.set_font('Arial', '', 12)
#        self.cell(0, 8, f'Page {self.page_no()}', 0, 0, 'C')
pdf = PDF()
pdf.add_page()

#line break
#pdf.ln(ch)

#add in the dRep primary dendrogram
#pdf.add_page()
#a description of the dendrogram
drepdescription="The primary clustering dendrogram output by dRep compare summarizes the pair-wise MASH distance between" \
                " all genomes in the genome list.\n The dotted line provides a visualization of the primary ANI threshold " \
                "(set to 80% for the unknownisolate automator). This tree is not created from genome alignments." \
                "This tree is created from the MASH-distances (distances between isolates calculated by MASH)" \
                "If you would like a phylogenetic tree, please run an " \
                "automator for tree-building (eg. bcgtree)."
pdf.set_font('Arial', 'B', 16)
pdf.cell(w=0, h=20, txt="Primary Clustering Dendrogram from dRep", ln=1)
pdf.image('Primary_clustering_dendrogram0.png',
          x = 10, y = None, w = 200, h = 200, type = 'PNG')
pdf.set_font('Arial', '', 12)
pdf.multi_cell(w=0, h=5, txt=drepdescription)

#add in some of the ANI matrices
#pdf_highlighter.process_data('{}'.format(args.ANIb_matrix),'ANIb_highlighted.pdf',"{}".format(args.seqid))

#description
anigraph_paragraph="Each cell represents a pariwise comparison between the named genomes on rows and columns. " \
                   "Cells with identity > 0.95 are red, those with <0.95 identity are blue. This division corresponds " \
                   "to a widely-used convention for bacterial species boundaries.\n " \
                   "(Information from: https://pyani.readthedocs.io/en/latest/interpreting_plots.html)"

#first add a page break
pdf.add_page()
pdf.set_font('Arial', 'B', 16)
pdf.cell(w=0, h=20, txt="Percentage identity matrix for ANIb analysis", ln=1)
pdf.set_font('Arial', '', 10)
pdf.multi_cell(w=0, h=5, txt="ANI (average nucleotide identity) is a pairwisemeasure of overall similarity between two genome sequences. \n"\
                             "ANIb analysis is considered the standard method in prokaryotic taxonomy. It uses BLASTN+ for average nucleotide identity " \
                             "calculations. The ANIb method aligns 1020nt fragments of the input sequences. \n"\
                             "More information regarding ANI calculations can be found at \n https://www.sciencedirect.com/science/article/pii/S0580951714000087")
pdf.image('{}'.format(args.ANIb_matrix),
#pdf.image('ANIb_highlighted.pdf', would need to change the automator to use the pdf file in order to highlight
          x = 10, y = None, w = 200, h = 200, type = 'PNG')
pdf.set_font('Arial', '', 10)
pdf.multi_cell(w=0, h=5, txt=anigraph_paragraph)

#do the same for ANIm
pdf.add_page()
pdf.set_font('Arial', 'B', 16)
pdf.cell(w=0, h=20, txt="Percentage identity matrix for ANIm analysis", ln=1)
pdf.set_font('Arial', '', 10)
pdf.multi_cell(w=0, h=5, txt="ANIm analysis uses MUMmer for average nucleotide identity calculations. The MUMmer program performs ultrafast "\
                             "alignment of two genomes. \n" \
                             "More information regarding ANI calculations can be found at \n https://www.sciencedirect.com/science/article/pii/S0580951714000087")
pdf.image('{}'.format(args.ANIm_matrix),
          x = 10, y = None, w = 200, h = 200, type = 'PNG')
pdf.set_font('Arial', '', 10)
pdf.multi_cell(w=0, h=5, txt=anigraph_paragraph)

#save output
pdf.output('GROBI_ANI_matrices.pdf', 'F')

#now use reportlab to create the fillable form coverpage
todaysdate= datetime.today().strftime('%Y-%m-%d')

def createSimpleForm():  
    my_canvas = canvas.Canvas('reportlab_coverpage.pdf', pagesize=letter)

    width, height = letter

    form = my_canvas.acroForm
      
    my_canvas.setFont("Helvetica-Bold", 20)  
    my_canvas.drawCentredString(300, 750, 'Genomic Report of Bacterial Identification (GROBI)')


    my_canvas.setFont("Helvetica", 16)  
    my_canvas.drawString(10, 700, 'Date: {}'.format(todaysdate))
    #contact email
    my_canvas.drawString(10, 670, 'Contact:')  
    form.textfield(  
        name = 'fname',  
        tooltip = 'Contact',  
        x = 70,  
        y = 665,  
        borderStyle = 'inset',  
        #borderColor = magenta,  
        #fillColor = pink,   
        width = 505,
        height = 20,  
        #textColor = blue,  
        forceBorder = True  
        )
    #my_canvas.drawString(10, 670, 'Contact: catherine.carrillo@inspection.gc.ca')
    

    my_canvas.setFont("Helvetica", 14)  
    #form = my_canvas.acroForm
      
    my_canvas.drawString(40, 640, 'Likely genus/species:')  
    form.textfield(  
        name = 'fname',  
        tooltip = 'Likely genus/species',  
        x = 175,  
        y = 625,  
        borderStyle = 'inset',  
        #borderColor = magenta,  
        #fillColor = pink,   
        width = 400,
        height = 30,  
        textColor = blue,  
        forceBorder = True  
        )

    my_canvas.setFont("Helvetica", 12)

    my_canvas.drawString(70, 600, 'Analyst signature:')  
    form.textfield(  
        name = 'fname',  
        tooltip = 'Analyst',  
        x = 175,  
        y = 590,  
        borderStyle = 'inset',  
        #borderColor = magenta,  
        #fillColor = pink,   
        width = 400,
        height = 25,  
        textColor = blue,  
        forceBorder = True  
        )

    my_canvas.drawString(132, 570, 'Notes:')  
    form.textfield(  
        name = 'fname',  
        tooltip = 'Notes',
        maxlen=1000,  
        x = 175,  
        y = 530,  
        borderStyle = 'inset',  
        #borderColor = magenta,  
        #fillColor = pink,   
        width = 400,  
        height = 50,
        #multiline=1,
        fieldFlags='multiline',
        textColor = blue,  
        forceBorder = True  
        )

    #add image table to first page
    grobiimage = 'GROBI_table.png'
    my_canvas.drawImage(grobiimage, 20, 100, width=600, preserveAspectRatio=True, mask='auto')

    my_canvas.showPage() #ends the page, sends things to next page

    my_canvas.save()

createSimpleForm()

#merge the coverpage and the ANI graphs into one pdf file
def merge_pdf_files_pdfrw(pdf_files, output_filename):
  output = PdfWriter()
  num = 0
  output_acroform = None
  for pdf in pdf_files:
      input = PdfReader(pdf,verbose=False)
      output.addpages(input.pages)
      if PdfName('AcroForm') in input[PdfName('Root')].keys():  # Not all PDFs have an AcroForm node
          source_acroform = input[PdfName('Root')][PdfName('AcroForm')]
          if PdfName('Fields') in source_acroform:
              output_formfields = source_acroform[PdfName('Fields')]
          else:
              output_formfields = []
          num2 = 0
          for form_field in output_formfields:
              key = PdfName('T')
              old_name = form_field[key].replace('(','').replace(')','')  # Field names are in the "(name)" format
              form_field[key] = 'FILE_{n}_FIELD_{m}_{on}'.format(n=num, m=num2, on=old_name)
              num2 += 1
          if output_acroform == None:
              # copy the first AcroForm node
              output_acroform = source_acroform
          else:
              for key in source_acroform.keys():
                  # Add new AcroForms keys if output_acroform already existing
                  if key not in output_acroform:
                      output_acroform[key] = source_acroform[key]
              # Add missing font entries in /DR node of source file
              if (PdfName('DR') in source_acroform.keys()) and (PdfName('Font') in source_acroform[PdfName('DR')].keys()):
                  if PdfName('Font') not in output_acroform[PdfName('DR')].keys():
                      # if output_acroform is missing entirely the /Font node under an existing /DR, simply add it
                      output_acroform[PdfName('DR')][PdfName('Font')] = source_acroform[PdfName('DR')][PdfName('Font')]
                  else:
                      # else add new fonts only
                      for font_key in source_acroform[PdfName('DR')][PdfName('Font')].keys():
                          if font_key not in output_acroform[PdfName('DR')][PdfName('Font')]:
                              output_acroform[PdfName('DR')][PdfName('Font')][font_key] = source_acroform[PdfName('DR')][PdfName('Font')][font_key]
          if PdfName('Fields') not in output_acroform:
              output_acroform[PdfName('Fields')] = output_formfields
          else:
              # Add new fields
              output_acroform[PdfName('Fields')] += output_formfields
      num +=1
  output.trailer[PdfName('Root')][PdfName('AcroForm')] = output_acroform
  output.write(output_filename)

merge_pdf_files_pdfrw(["reportlab_coverpage.pdf","GROBI_ANI_matrices.pdf"],'{}'.format(args.outfile))

