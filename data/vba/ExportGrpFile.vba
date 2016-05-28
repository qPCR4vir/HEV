Attribute VB_Name = "ExportGrpFile"
Sub SelectThisSeqRegion()
Attribute SelectThisSeqRegion.VB_ProcData.VB_Invoke_Func = "r\n14"
'
' SelectThisSeqRegion Makro
'
' Tastenkombination: Strg+r
'
' Let you want to analyze a sequence that you have just enter in the big alignment.
' You have already inserted a row in the large table with the correct name (the same as in the alignment used in MEGA)
' and with the position of biginning and end of the sequence in the alignment.
' You now need to set the region and export the column "span" into a group file.
' With this macro you just click on the sequences row to select it entirely and then press Ctrl-r.
' This macro first transfers the beginning and end of the selected sequence to the sheet "Region",
' and sets the option "manual selection" there to trigger the selection of all the sequences
' which span the same region of the alignment.
' Then it copy the "span" column to a new workbook for which may create a directory Aut.MEGA\
' if not yet exist.
' Before closing the new workbook, it is saved as a text file with a name that show the region selected.
    
    Application.CutCopyMode = False
    Dim beg As Integer, fin As Integer, FileName As String
     
    beg = Application.Intersect(Selection, Range("SeqBeg")).Value ' where the seq begin in the aligment
    fin = Application.Intersect(Selection, Range("SeqEnd")).Value ' where the seq end
    Range("manual_beg").Value = beg    ' set as the manually selected region
    Range("manual_end").Value = fin
    Range("GenRegSelect").Value = 3    ' select the manually selected region in the big table
                                       ' this generate the column with the "span" or "out"
    Range("span_column").Copy
    
    
    FileName = Application.ActiveWorkbook.Path & "\HEV." & Str(beg) & "-" & Str(fin) & ".grp.txt"

    With Workbooks.Add   ' create a new workbook file
      .Sheets(1).Range("A1").PasteSpecial Paste:=xlPasteValues, Operation:=xlNone, SkipBlanks:=False, Transpose:=False
      .SaveAs FileName:=FileName, FileFormat:=xlTextMSDOS
      .Close SaveChanges:=False
    End With
    
    MsgBox ("A file " & FileName & " have been created.")
      
End Sub

