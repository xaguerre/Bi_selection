#ifndef SNDISPLAY_CC
#define SNDISPLAY_CC

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TH2D.h"
#include "TLine.h"
#include "TPaletteAxis.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TString.h"
#include "TText.h"

#include<vector>

namespace sndisplay
{
  ////////////////////////
  // sndisplay::palette //
  ////////////////////////

  class palette
  {
  private:
    palette()
    {
      const Int_t nRGBs = 6;
      Double_t stops[nRGBs] = { 0.00, 0.20, 0.40, 0.60, 0.80, 1.00 };
      Double_t red[nRGBs]   = { 0.25, 0.00, 0.20, 1.00, 1.00, 0.90 };
      Double_t green[nRGBs] = { 0.25, 0.80, 1.00, 1.00, 0.80, 0.00 };
      Double_t blue[nRGBs]  = { 1.00, 1.00, 0.20, 0.00, 0.00, 0.00 };

      palette_index = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 100);
    }

    static palette *instance;
    int palette_index;

  public:
    ~palette() {};

    static palette *get_me() {
      if (instance == nullptr) instance = new palette;
      return instance;}

    static int get_index() {
      return get_me()->palette_index;}

  }; // sndisplay::palette class

  palette *palette::instance = nullptr;

  ///////////////////////
  // sndisplay::canvas //
  ///////////////////////

  // wrapper to create a TCanvas with exact width/height
  // (taking into account OS window decorator) and with
  // a fixed aspect ratio

  TCanvas* canvas (const char *name, const char *title, int width, int height)
  {
    TCanvas *c = new TCanvas (name, title, width, height);

    int decoration_width = width - c->GetWw();
    int decoration_height = height - c->GetWh();

    c->SetWindowSize(width+decoration_width, height+decoration_height);

    c->SetFixedAspectRatio();

    return c;
  }

  ////////////////////////////
  // sndisplay::calorimeter //
  ////////////////////////////

  class calorimeter
  {
  public:
    calorimeter (const char *n = "", bool with_p=false) : calorimeter_name (n), with_palette(with_p)
    {
      // 1 canvas mode
      canvas = nullptr;

      // 2 canvas mode
      canvas_it = nullptr;
      canvas_fr = nullptr;

      pad_it = nullptr;
      pad_fr = nullptr;
      pad_palette = nullptr;
      pad_palette_it = nullptr;
      pad_palette_fr = nullptr;

      // draw by default the OM ID
      draw_omid = true;
      draw_omnum = false;
      draw_content = false;
      draw_content_format = "%.0f";

      for (int omnum=0; omnum<nb_om; ++omnum)
	content.push_back(0);

      range_min = range_max = -1;

      const double spacerx = 0.01; // with_palette ? 0.0093458 : 0.0100;
      const double spacery = 0.0125;

      const double mw_sizey = (1-4*spacery)/(13+2);
      const double gv_sizey = mw_sizey;
      const double xw_sizey = mw_sizey*13/16.;

      const double mw_sizex = (1-4*spacerx)/(20+4); // with_palette ? (1-5*spacerx)/(20+4+1.5) : (1-4*spacerx)/(20+4);
      const double gv_sizex = mw_sizex*20./16.;
      const double xw_sizex = mw_sizex;

      //////////////////////////
      // MWALL initialisation //
      //////////////////////////

      for (int mw_side=0; mw_side<2; ++mw_side) {

	for (int mw_column=0; mw_column<20; ++mw_column) {

	  for (int mw_row=0; mw_row<13; ++mw_row) {

	    int omnum = mw_side*20*13 + mw_column*13 + mw_row;

	    double x1 = spacerx + 2*xw_sizex + spacerx; // + 0.5*mw_side;
	    x1 += (mw_side == 0) ? mw_sizex*(19-mw_column) : mw_sizex*(mw_column); // swap IT for external view

	    double y1 = spacery + gv_sizey + spacery + mw_sizey*(mw_row);
	    double x2 = x1 + mw_sizex;
	    double y2 = y1 + mw_sizey;

	    TBox box (x1, y1, x2, y2);
	    box.SetFillColor(0);
	    box.SetLineWidth(1);
	    ombox.push_back(box);

	    TString omid_string = Form("M:%1d.%d.%d", mw_side, mw_column, mw_row);
	    TText omid_text (x1+0.5*mw_sizex, y1+0.7*mw_sizey, omid_string);
	    omid_text.SetTextFont(42);
	    omid_text.SetTextSize(0.013);
	    omid_text.SetTextAlign(22);
	    omid_text_v.push_back(omid_text);

	    TString omnum_string = Form("%03d", omnum);
	    TText omnum_text (x1+0.5*mw_sizex, y1+0.7*mw_sizey, omnum_string);
	    omnum_text.SetTextFont(42);
	    omnum_text.SetTextSize(0.013);
	    omnum_text.SetTextAlign(22);
	    omnum_text_v.push_back(omnum_text);

	    TText content_text (x1+0.5*mw_sizex, y1+0.3*mw_sizey, "");
	    content_text.SetTextSize(0.02);
	    content_text.SetTextAlign(22);
	    content_text_v.push_back(content_text);

	  } // for mw_row

	} // for mw_column

      } // for mw_side

      //////////////////////////
      // XWALL initialisation //
      //////////////////////////

      for (int xw_side=0; xw_side<2; ++xw_side) {

	for (int xw_wall=0; xw_wall<2; ++xw_wall) {

	  for (int xw_column=0; xw_column<2; ++xw_column) {

	    for (int xw_row=0; xw_row<16; ++xw_row) {

	    int omnum = 520 + xw_side*2*2*16 + xw_wall*2*16 + xw_column*16 + xw_row;

	    double x1;

	    switch (xw_side) {
	    case 0:
	      if (xw_wall == 0) // IT-MO
		x1 = spacerx + 2*xw_sizex + spacerx + 20*mw_sizex + spacerx + (1-xw_column)*xw_sizex;
	      else              // IT-TU
		x1 = spacerx + xw_sizex*xw_column;
	      break;

	    case 1:
	      if (xw_wall == 0) // FR-MO
		// x1 = 0.5 + spacerx + xw_sizex*xw_column;
		x1 = spacerx + xw_sizex*xw_column;
	      else // FR-TU
		// x1 = 0.5 + spacerx + 2*xw_sizex + spacerx + 20*mw_sizex + spacerx + (1-xw_column)*xw_sizex;
		x1 = spacerx + 2*xw_sizex + spacerx + 20*mw_sizex + spacerx + (1-xw_column)*xw_sizex;
	      break;}

	    double x2 = x1 + xw_sizex;
	    double y1 = spacery + gv_sizey + spacery + xw_sizey*(xw_row);
	    double y2 = spacery + gv_sizey + spacery + xw_sizey*(xw_row+1);

	    TBox box (x1, y1, x2, y2);
	    box.SetFillColor(0);
	    box.SetLineWidth(1);
	    ombox.push_back(box);

	    TString omid_string = Form("X:%1d.%1d.%1d.%d", xw_side, xw_wall, xw_column, xw_row);

	    TText omid_text (x1+0.5*mw_sizex, y1+0.7*xw_sizey, omid_string);
	    omid_text.SetTextFont(42);
	    omid_text.SetTextSize(0.013);
	    omid_text.SetTextAlign(22);
	    omid_text_v.push_back(omid_text);

	    TString omnum_string = Form("%03d", omnum);

	    TText omnum_text (x1+0.5*xw_sizex, y1+0.7*xw_sizey, omnum_string);
	    omnum_text.SetTextFont(42);
	    omnum_text.SetTextSize(0.013);
	    omnum_text.SetTextAlign(22);
	    omnum_text_v.push_back(omnum_text);

	    TText content_text (x1+0.5*xw_sizex, y1+0.3*xw_sizey, "");
	    content_text.SetTextSize(0.02);
	    content_text.SetTextAlign(22);
	    content_text_v.push_back(content_text);

	    } // for xw_row

	  } // for xw_column

	} // for xw_wall

      } // for xw_side

      //////////////////////////
      // GVETO initialisation //
      //////////////////////////

      for (int gv_side=0; gv_side<2; ++gv_side) {

	for (int gv_wall=0; gv_wall<2; ++gv_wall) {

	  for (int gv_column=0; gv_column<16; ++gv_column) {

	    int omnum = 520 + 128 + gv_side*2*16 + gv_wall*16 + gv_column;

	    double x1;

	    if (gv_side == 0)
	      x1 = spacerx + 2*xw_sizex + spacerx + gv_sizex*(16-1-gv_column);
	    else
	      // x1 = 0.5 + spacerx + 2*xw_sizex + spacerx + gv_sizex*gv_column;
	      x1 = spacerx + 2*xw_sizex + spacerx + gv_sizex*gv_column;

	    double x2 = x1 + gv_sizex;
	    double y1 = spacery + gv_wall*(gv_sizey + spacery + 13*mw_sizey + spacery);
	    double y2 = y1 + gv_sizey;

	    TBox box (x1, y1, x2, y2);
	    box.SetFillColor(0);
	    box.SetLineWidth(1);
	    ombox.push_back(box);

	    TString omid_string = Form("G:%1d.%1d.%d", gv_side, gv_wall, gv_column);

	    TText omid_text (x1+0.5*gv_sizex, y1+0.7*gv_sizey, omid_string);
	    omid_text.SetTextFont(42);
	    omid_text.SetTextSize(0.013);
	    omid_text.SetTextAlign(22);
	    omid_text_v.push_back(omid_text);

	    TString omnum_string = Form("%03d", omnum);

	    TText omnum_text (x1+0.5*gv_sizex, y1+0.7*gv_sizey, omnum_string);
	    omnum_text.SetTextFont(42);
	    omnum_text.SetTextSize(0.013);
	    omnum_text.SetTextAlign(22);
	    omnum_text_v.push_back(omnum_text);

	    TText content_text (x1+0.5*gv_sizex, y1+0.3*gv_sizey, "");
	    content_text.SetTextSize(0.02);
	    content_text.SetTextAlign(22);
	    content_text_v.push_back(content_text);

	  } // for gv_column

	} // for gv_wall

      } // for gv_side

      label_it = new TText (spacerx+xw_sizex, spacery+gv_sizey+spacery+13*mw_sizey+spacery+0.5*gv_sizey, "  ITALY");
      label_it->SetTextSize(0.028);
      label_it->SetTextAlign(22);

      label_fr = new TText (spacerx+xw_sizex, spacery+gv_sizey+spacery+13*mw_sizey+spacery+0.5*gv_sizey, "FRANCE");
      label_fr->SetTextSize(0.028);
      label_fr->SetTextAlign(22);

      palette_histo = nullptr;
      palette_axis = nullptr;

      const double palette_sizey = mw_sizey*13;
      // const double palette_sizex = mw_sizex*1.5;

      palette_histo = new TH2D(Form("%s_palette_histo",calorimeter_name.Data()), "", 1, 0, 1, 1, 0, 1);
      palette_histo->GetZaxis()->SetLabelSize(0.24);
      palette_histo->GetZaxis()->SetLabelOffset(0.1);
      palette_histo->GetZaxis()->SetNdivisions(509);
      palette_histo->GetZaxis()->SetTickLength(0.65);
      // palette_histo->GetZaxis()->SetLabelFont(62);
      palette_histo->SetMinimum(range_min);
      palette_histo->SetMaximum(range_max);
      palette_histo->SetContour(100);

      // the constructor TPaletteAxis(x1, y1, x2, y2, histo)
      // crashing due to no canvas existing (gPad = nullptr)
      palette_axis = new TPaletteAxis;
      palette_axis->SetHistogram(palette_histo);
      palette_axis->SetX1NDC(0.1); // 1-spacerx-palette_sizex*(5./6)); // position tuning
      palette_axis->SetY1NDC(spacery+gv_sizey+spacery+palette_sizey/8.);
      palette_axis->SetX2NDC(0.5); //1-spacerx-palette_sizex*(3./6)); // position tuning
      palette_axis->SetY2NDC(1-palette_axis->GetY1NDC());

    }; // calorimeter()

    ~calorimeter()
    {
      // if (palette_axis) delete palette_axis;
      // if (palette_histo) delete palette_histo;

      delete label_fr;
      delete label_it;

      delete canvas_fr;
      delete canvas_it;
    };

    static const int nmwall = 520;
    static const int nxwall = 128;
    static const int ngveto =  64;

    static const int nb_om  = 712;

    void setrange(float zmin, float zmax)
    {
      range_min = zmin; range_max = zmax;
    }

    void draw_omid_label (bool draw=true)
    {
      draw_omid = draw;
    }

    void draw_omnum_label(bool draw=true)
    {
      draw_omnum = draw;
    }

    void draw_content_label(const char *format="%.0f")
    {
      draw_content_format = TString(format);
      draw_content = true;
    }

    // draw in 1 canvas mode
    void draw1()
    {
      if (pad_palette_fr != nullptr) {delete pad_palette_fr; pad_palette_fr = nullptr;}
      if (pad_palette_it != nullptr) {delete pad_palette_it; pad_palette_it = nullptr;}
      if (pad_palette != nullptr) {delete pad_palette; pad_palette = nullptr;}
      if (pad_fr != nullptr) {delete pad_fr; pad_it = nullptr;}
      if (pad_it != nullptr) {delete pad_it; pad_it = nullptr;}
      if (canvas_fr != nullptr) {delete canvas_fr; canvas_fr = nullptr;}
      if (canvas_it != nullptr) {delete canvas_it; canvas_it = nullptr;}
      if (canvas != nullptr) {delete canvas; canvas = nullptr;}

      const int calorimeter_width = 1200;
      const int palette_width = 42;

      const int canvas_width  = with_palette ? (calorimeter_width + palette_width) : calorimeter_width;
      const int canvas_height = 390;

      canvas = sndisplay::canvas(Form("%s_canvas",calorimeter_name.Data()), Form("%s",calorimeter_name.Data()), canvas_width, canvas_height);

      if (with_palette)
	{
	  const double palette_rwidth = ((double)(palette_width))/((double)(canvas_width));
	  const double pad_x1 = 0.5 - palette_rwidth/2;
	  const double pad_x2 = 1.0 - palette_rwidth;

	  pad_it = new TPad (Form("%s_pad_it",calorimeter_name.Data()), Form("%s it",calorimeter_name.Data()), 0, 0, pad_x1, 1);
	  pad_fr = new TPad (Form("%s_pad_fr",calorimeter_name.Data()), Form("%s fr",calorimeter_name.Data()), pad_x1, 0, pad_x2, 1);
	  pad_palette = new TPad (Form("%s_pad_palette",calorimeter_name.Data()), Form("%s palette",calorimeter_name.Data()), pad_x2, 0, 1, 1);
	}
      else
	{
	  pad_it = new TPad (Form("%s_pad_it",calorimeter_name.Data()), Form("%s it",calorimeter_name.Data()), 0, 0, 0.5, 1);
	  pad_fr = new TPad (Form("%s_pad_fr",calorimeter_name.Data()), Form("%s fr",calorimeter_name.Data()), 0.5, 0, 1, 1);
	}

      canvas->cd();
      pad_it->Draw();
      pad_fr->Draw();

      draw_internal_content();

      if (with_palette)
	{
	  canvas->cd();
	  pad_palette->Draw();
	  pad_palette->cd();
	  palette_axis->Draw();
	}

      canvas->SetEditable(false);
    }

    // draw in 2 canvas mode
    void draw2()
    {
      if (pad_palette_fr != nullptr) {delete pad_palette_fr; pad_palette_fr = nullptr;}
      if (pad_palette_it != nullptr) {delete pad_palette_it; pad_palette_it = nullptr;}
      if (pad_palette != nullptr) {delete pad_palette; pad_palette = nullptr;}
      if (pad_fr != nullptr) {delete pad_fr; pad_it = nullptr;}
      if (pad_it != nullptr) {delete pad_it; pad_it = nullptr;}
      if (canvas_fr != nullptr) {delete canvas_fr; canvas_fr = nullptr;}
      if (canvas_it != nullptr) {delete canvas_it; canvas_it = nullptr;}
      if (canvas != nullptr) {delete canvas; canvas = nullptr;}

      const int calorimeter_width = 1200;
      const int palette_width = 84;

      const int canvas_width  = with_palette ? (calorimeter_width + palette_width) : calorimeter_width;
      const int canvas_height = 780;

      const double pad_xmax = with_palette ? ((double)(calorimeter_width))/((double)(canvas_width)) : 1;

      canvas_it = sndisplay::canvas (Form("%s_canvas_it",calorimeter_name.Data()), Form("%s (IT side)",calorimeter_name.Data()), canvas_width, canvas_height);
      pad_it = new TPad (Form("%s_pad_it",calorimeter_name.Data()), Form("%s it",calorimeter_name.Data()), 0, 0, pad_xmax, 1);
      pad_it->Draw();

      canvas_fr = sndisplay::canvas (Form("%s_canvas_fr",calorimeter_name.Data()), Form("%s (FR side)",calorimeter_name.Data()), canvas_width, canvas_height);
      pad_fr = new TPad (Form("%s_pad_fr",calorimeter_name.Data()), Form("%s fr",calorimeter_name.Data()), 0, 0, pad_xmax, 1);
      pad_fr->Draw();

      if (with_palette)
	{
	  pad_palette_it = new TPad (Form("%s_pad_palette_it",calorimeter_name.Data()), Form("%s palette it",calorimeter_name.Data()), pad_xmax, 0, 1, 1);
	  pad_palette_fr = new TPad (Form("%s_pad_palette_fr",calorimeter_name.Data()), Form("%s palette fr",calorimeter_name.Data()), pad_xmax, 0, 1, 1);
	}

      draw_internal_content();

      if (with_palette)
	{
	  canvas_it->cd();
	  pad_palette_it->Draw();
	  pad_palette_it->cd();
	  palette_axis->Draw();

	  canvas_fr->cd();
	  pad_palette_fr->Draw();
	  pad_palette_fr->cd();
	  palette_axis->Draw();
	}

      canvas_it->SetEditable(false);
      canvas_fr->SetEditable(false);
    }

    void draw_internal_content()
    {
      // update color and content
      update(false);

      /////////////
      // Draw IT //
      /////////////

      pad_it->cd();

      int mw_side=0;
      for (int mw_column=0; mw_column<20; ++mw_column) {
	for (int mw_row=0; mw_row<13; ++mw_row) {
	  int id = mw_side*20*13 + mw_column*13 + mw_row;
	  ombox[id].Draw("l");
	  if (draw_omid)
	    omid_text_v[id].Draw();
	  else if (draw_omnum)
	    omnum_text_v[id].Draw();
	  if ((draw_content && content[id]!=0) || text_was_set)
	    content_text_v[id].Draw();
	}
      }

      int xw_side=0;
      for (int xw_wall=0; xw_wall<2; ++xw_wall) {
	for (int xw_column=0; xw_column<2; ++xw_column) {
	  for (int xw_row=0; xw_row<16; ++xw_row) {
	    int id = 520 + xw_side*2*2*16 + xw_wall*2*16 + xw_column*16 + xw_row;
	    ombox[id].Draw("l");
	    if (draw_omid)
	      omid_text_v[id].Draw();
	    else if (draw_omnum)
	      omnum_text_v[id].Draw();
	    if ((draw_content && content[id]!=0) || text_was_set)
	      content_text_v[id].Draw();
	  }
	}
      }

      int gv_side=0;
      for (int gv_wall=0; gv_wall<2; ++gv_wall) {
	for (int gv_column=0; gv_column<16; ++gv_column) {
	  int id = 520 + 128 + gv_side*2*16 + gv_wall*16 + gv_column;
	  ombox[id].Draw("l");
	  if (draw_omid) // if ((gv_column % 5) == 0)
	    omid_text_v[id].Draw();
	  else if (draw_omnum)
	    omnum_text_v[id].Draw();
	  if ((draw_content && content[id]!=0) || text_was_set)
	    content_text_v[id].Draw();
	}
      }

      label_it->Draw();

      /////////////
      // Draw FR //
      /////////////

      pad_fr->cd();

      mw_side=1;
      for (int mw_column=0; mw_column<20; ++mw_column) {
	for (int mw_row=0; mw_row<13; ++mw_row) {
	  int id = mw_side*20*13 + mw_column*13 + mw_row;
	  ombox[id].Draw("l");
	  if (draw_omid)
	      omid_text_v[id].Draw();
	  else if (draw_omnum)
	    omnum_text_v[id].Draw();
	  if (draw_content || text_was_set)
	    content_text_v[id].Draw();
	}
      }

      xw_side=1;
      for (int xw_wall=0; xw_wall<2; ++xw_wall) {
	for (int xw_column=0; xw_column<2; ++xw_column) {
	  for (int xw_row=0; xw_row<16; ++xw_row) {
	    int id = 520 + xw_side*2*2*16 + xw_wall*2*16 + xw_column*16 + xw_row;
	    ombox[id].Draw("l");
	    if (draw_omid)
	      omid_text_v[id].Draw();
	    else if (draw_omnum)
	      omnum_text_v[id].Draw();
	    if (draw_content || text_was_set)
	      content_text_v[id].Draw();
	  }
	}
      }

      gv_side=1;
      for (int gv_wall=0; gv_wall<2; ++gv_wall) {
	for (int gv_column=0; gv_column<16; ++gv_column) {
	  int id = 520 + 128 + gv_side*2*16 + gv_wall*16 + gv_column;
	  ombox[id].Draw("l");
	  if (draw_omid)
	    omid_text_v[id].Draw();
	  else if (draw_omnum)
	    omnum_text_v[id].Draw();
	  if (draw_content || text_was_set)
	    content_text_v[id].Draw();
	}
      }

      label_fr->Draw();

      text_was_set = false;

      update_canvas();
    }

    void draw()
    {
      draw2();
    }

    void palette (bool with_p=true)
    {
      if (with_palette == with_p)
	return; // do nothing

      with_palette = with_p;

      if (canvas != nullptr) draw1();
      else if (canvas_it != nullptr) draw2();
      else draw1();
    }

    void reset()
    {
      for (int omnum=0; omnum<nb_om; ++omnum)
	{
	  content[omnum] = 0;
	  content_text_v[omnum].Clear();
	  ombox[omnum].SetFillColor(0);
	}
    }

    float getcontent (int omnum)
    {
      return content[omnum];
    }

    void setcontent (int omnum, float value)
    {
      if ((omnum >= 0) && (omnum < nb_om))
	content[omnum] = value;
      else printf("*** wrong OM NUM\n");
    }

    void setcontent (int om_side, int om_wall, int om_column, int om_row, float value)
    {
      int omnum = -1;

      // auto detect MW
      if ((om_side!=-1) && (om_wall==-1) && (om_column!=-1) && (om_row!=-1))
	omnum = 260*om_side + 13*om_column + om_row;

      // auto detect XW
      else if ((om_side!=-1) && (om_wall!=-1) && (om_column!=-1) && (om_row!=-1))
	omnum = 520 + 64*om_side + 32*om_wall + 16*om_column + om_row;

      // auto detect GV
     else if ((om_side!=-1) && (om_wall!=-1) && (om_column!=-1) && (om_row==-1))
	omnum = 520 + 128 + 32*om_side + 16*om_wall + om_column;

      else {
	printf("+++ sndisplay: skipping OM (%d.%d.%d.%d)\n", om_side, om_wall, om_column, om_row);
	return;}

      content[omnum] = value;
    }

    void setcolor (int omnum, Color_t color)
    {
      if ((omnum >= 0) && (omnum < nb_om))
	{
	  ombox[omnum].SetFillColor(color);
	  color_was_set = true;
	}
      else printf("*** wrong OM NUM\n");
    }

    void settext (int omnum, const char *text)
    {
      if ((omnum >= 0) && (omnum < nb_om))
	{
	  content_text_v[omnum].SetText(content_text_v[omnum].GetX(), content_text_v[omnum].GetY(), text);
	  text_was_set = true;
	}
      else printf("*** wrong OM NUM\n");
    }

    void setmwcolor (int om_side, int om_column, int om_row, Color_t color)
    {
      setcolor(260*om_side + 13*om_column + om_row, color);
    }

    void setmwtext (int om_side, int om_column, int om_row, const char *text)
    {
      settext(260*om_side + 13*om_column + om_row, text);
    }

    void setxwcolor (int om_side, int om_wall, int om_column, int om_row, Color_t color)
    {
      setcolor(520 + 64*om_side + 32*om_wall + 16*om_column + om_row, color);
    }

    void setxwtext (int om_side, int om_wall, int om_column, int om_row, const char *text)
    {
      settext(520 + 64*om_side + 32*om_wall + 16*om_column + om_row, text);
    }

    void setgvcolor (int om_side, int om_wall, int om_column, Color_t color)
    {
      setcolor(648 + 32*om_side + 16*om_wall + om_column, color);
    }

    void setgvtext (int om_side, int om_wall, int om_column, const char *text)
    {
      settext(648 + 32*om_side + 16*om_wall + om_column, text);
    }

    void fill (int omnum, float value=1)
    {
      setcontent(omnum, content[omnum]+value);
    }

    void update_canvas()
    {
      if (canvas != nullptr)
	{
	  canvas->Modified();
	  canvas->Update();
	}

      if (canvas_it != nullptr)
	{
	  canvas_it->Modified();
	  canvas_it->Update();
	}

      if (canvas_fr != nullptr)
	{
	  canvas_fr->Modified();
	  canvas_fr->Update();
	}

      gSystem->ProcessEvents();
    }

    void update (bool update_canvas_too=true)
    {
      // skip the automated color range by OM content
      // if some OM color was set with setxxcolor()

      if (color_was_set)
	{
	  color_was_set = false;
	  return;
	}

      // autoset Z range [0, content_max]
      // unless setrange() has been called

      float content_min = content[0];
      float content_max = content[0];

      for (int omnum=1; omnum<nb_om; ++omnum)
	{
	  if (content[omnum] < content_min) content_min = content[omnum];
	  if (content[omnum] > content_max) content_max = content[omnum];
	}

      content_min = 0;

      if (range_min != -1) content_min = range_min;
      if (range_max != -1) content_max = range_max;

      // update color box

      for (int omnum=0; omnum<nb_om; ++omnum)
	{
	  if (content[omnum] != 0)
	    {
	      int color_index = floor (99*(content[omnum]-content_min)/(content_max-content_min));
	      if (color_index < 0) color_index = 0;
	      else if (color_index >= 100) color_index = 99;
	      ombox[omnum].SetFillColor(palette::get_index() + color_index);
	    }
	  else
	    ombox[omnum].SetFillColor(0);
	}

      // update content text

      if (draw_content && !text_was_set)
	{
	  for (int omnum=0; omnum<nb_om; ++omnum)
	    {
	      TString text = "";
	      if (content[omnum] != 0)
		text = Form(draw_content_format.Data(), content[omnum]);
	      content_text_v[omnum].SetText(content_text_v[omnum].GetX(), content_text_v[omnum].GetY(), text);
	    }
	}

      // update palette

      if (palette_axis)
	{
	  palette_histo->SetMinimum(content_min);
	  palette_histo->SetMaximum(content_max);
	  palette_histo->SetContour(100);
	}

      //

      if (update_canvas_too)
	update_canvas();
    }

    TString calorimeter_name;

    // draw options
    bool with_palette;
    bool draw_omid;
    bool draw_omnum;
    bool draw_content;
    TString draw_content_format;
    float range_min, range_max;

    bool color_was_set;
    bool text_was_set;

    TCanvas *canvas;
    TCanvas *canvas_it;
    TCanvas *canvas_fr;

    TPad *pad_it;
    TPad *pad_fr;
    TPad *pad_palette;
    TPad *pad_palette_it;
    TPad *pad_palette_fr;

    TText *label_it;
    TText *label_fr;

    TH2D *palette_histo;
    TPaletteAxis *palette_axis;

    std::vector<float> content;

    std::vector<TBox>  ombox;
    std::vector<TText> omid_text_v;
    std::vector<TText> omnum_text_v;
    std::vector<TText> content_text_v;

  }; // sndisplay::calorimeter class

  ////////////////////////
  // sndisplay::tracker //
  ////////////////////////

  class tracker
  {
  public:
    tracker (const char *n = "", bool with_palette=false) : tracker_name (n)
    {
      canvas = nullptr;

      // do not draw by default the CELL ID
      draw_cellid = false;
      draw_cellnum = false;
      draw_content = false;
      draw_content_format = "%.0f";

      for (int cellnum=0; cellnum<nb_cell; ++cellnum)
	content.push_back(0);

      range_min = range_max = -1;

      const double spacerx = with_palette ? 0.0082 : 0.0082;
      const double spacery = 0.0500;

      const double cell_sizex = with_palette ? (1-3*spacerx)/(113.+2) : (1-2*spacerx)/113.0;
      const double cell_sizey = (1-3*spacery)/(9*2);

      //////////////////////////
      // CELLS initialisation //
      //////////////////////////

      for (int cell_side=0; cell_side<2; ++cell_side) {

	  for (int cell_row=0; cell_row<113; ++cell_row) {

	    double x1 = spacerx + cell_row*cell_sizex;

	    for (int cell_layer=0; cell_layer<9; ++cell_layer) {

	    int cellnum = cell_side*113*9 + cell_row*9 + cell_layer;

	    double y1 = spacery;

	    if (cell_side == 0)
	      y1 += 9*cell_sizey + spacery + cell_layer*cell_sizey;
	    else y1 += (8-cell_layer)*cell_sizey;

	    double x2 = x1 + cell_sizex;
	    double y2 = y1 + cell_sizey;

	    TBox box (x1, y1, x2, y2);
	    box.SetFillColor(0);
	    box.SetLineWidth(1);
	    cellbox.push_back(box);

	    TString cellid_string = Form("M:%1d.%d.%d", cell_side, cell_layer, cell_row);
	    TText cellid_text (x1+0.5*cell_sizex, y1+0.667*cell_sizey, cellid_string);
	    cellid_text.SetTextSize(0.01);
	    cellid_text.SetTextAlign(22);
	    cellid_text_v.push_back(cellid_text);

	    TString cellnum_string = Form("%03d", cellnum);
	    TText cellnum_text (x1+0.5*cell_sizex, y1+0.333*cell_sizey, cellnum_string);
	    cellnum_text.SetTextSize(0.01);
	    cellnum_text.SetTextAlign(22);
	    cellnum_text_v.push_back(cellnum_text);

	    TText content_text (x1+0.5*cell_sizex, y1+0.333*cell_sizey, "");
	    content_text.SetTextSize(0.01);
	    content_text.SetTextAlign(22);
	    content_text_v.push_back(content_text);

	    } // for cell_layer

	    if ((cell_row %5) == 0)
	      {
		TText row_text (x1+0.5*cell_sizex, 0.5, Form("%d",cell_row));
		row_text.SetTextSize(0.03);
		row_text.SetTextAngle(90);
		row_text.SetTextAlign(22);
		row_text_v.push_back(row_text);
	      }
	  } // for cell_row

      } // for cell_side

      label_it = new TText (spacerx, 2.5*spacery+2*9*cell_sizey, "ITALY");
      label_it->SetTextSize(0.036);
      label_it->SetTextAlign(12);

      // label_fr = new TText (0.5 + spacerx, spacery+gv_sizey+spacery+13*mw_sizey+spacery+0.25*gv_sizey, "FRANCE");
      label_fr = new TText (spacerx, 0.5*spacery, "FRANCE");
      label_fr->SetTextSize(0.036);
      label_fr->SetTextAlign(12);

      palette_histo = nullptr;
      palette_axis = nullptr;

      if (with_palette)
	{
	  const double palette_sizex = 2*cell_sizex;
	  const double palette_sizey = cell_sizey*9*2 + spacery;

	  palette_histo = new TH2D(Form("%s_palette_histo",tracker_name.Data()), "", 1, 0, 1, 1, 0, 1);
	  palette_histo->GetZaxis()->SetNdivisions(509);
	  palette_histo->GetZaxis()->SetLabelSize(0.032);
	  palette_histo->GetZaxis()->SetLabelFont(62);
	  palette_histo->GetZaxis()->SetTickLength(0.009);
	  palette_histo->SetMinimum(range_min);
	  palette_histo->SetMaximum(range_max);
	  palette_histo->SetContour(100);

	  // the constructor TPaletteAxis(x1, y1, x2, y2, histo)
	  // crashing due to no canvas existing (gPad = nullptr)
	  palette_axis = new TPaletteAxis;
	  palette_axis->SetHistogram(palette_histo);
	  palette_axis->SetX1NDC(1-spacerx-palette_sizex); // position tuning
	  palette_axis->SetY1NDC(spacery);
	  palette_axis->SetX2NDC(1-spacerx-palette_sizex/2); // position tuning
	  palette_axis->SetY2NDC(spacery+palette_sizey);
	  palette_axis->SetY2NDC(1-palette_axis->GetY1NDC());
	}

    }; // tracker()

    ~tracker()
    {
      // if (palette_axis) delete palette_axis;
      // if (palette_histo) delete palette_histo;

      delete label_fr;
      delete label_it;

      delete canvas;
    }

    static const int nb_cell  = 2034;

    void setrange(float zmin, float zmax)
    {
      range_min = zmin; range_max = zmax;
    }

    void draw_cellid_label() {
      draw_cellid = true;}

    void draw_cellnum_label() {
      draw_cellnum = true;}

    void draw_content_label(const char *format="%.0f") {
      draw_content_format = TString(format);
      draw_content = true;}

    void draw()
    {
      const int canvas_width  = palette_axis ? 1231*2 : 1200*2;
      const int canvas_height =  221*2;

      update(false);

      if (canvas == nullptr)
	{
	  canvas = new TCanvas (Form("%s_canvas",tracker_name.Data()), tracker_name, canvas_width/2, canvas_height/2);

	  // force canvas exact size
	  int decoration_width = canvas_width/2 - canvas->GetWw();
	  int decoration_height = canvas_height/2 - canvas->GetWh();
	  int scroll_height = 16; // found by hand !

	  // adjust the canvas width to cover 3/4 of tracker width
	  canvas->SetWindowSize(canvas_width*3/4+decoration_width, canvas_height+decoration_height+scroll_height);
	  canvas->SetCanvasSize(canvas_width, canvas_height);
	}

      if (draw_content)
	{
	  for (int cellnum=0; cellnum<nb_cell; ++cellnum)
	    content_text_v[cellnum].SetText(content_text_v[cellnum].GetX(), content_text_v[cellnum].GetY(), Form(draw_content_format.Data(), content[cellnum]));
	}

      canvas->cd();
      canvas->SetEditable(true);

      for (int cell_side=0; cell_side<2; ++cell_side) {

	for (int cell_row=0; cell_row<113; ++cell_row) {

	  for (int cell_layer=0; cell_layer<9; ++cell_layer) {

	    int cellnum = cell_side*113*9 + cell_row*9 + cell_layer;

	    cellbox[cellnum].Draw("l");

	    if (draw_cellid)
	      cellid_text_v[cellnum].Draw();

	    else if (draw_cellnum)
	      cellnum_text_v[cellnum].Draw();

	    if ((draw_content && content[cellnum]!=0) || text_was_set)
	      content_text_v[cellnum].Draw();
	  }

	}

      }

      if (palette_axis) palette_axis->Draw();

      for (size_t row = 0; row<row_text_v.size(); ++row)
	row_text_v[row].Draw();

      label_it->Draw();
      label_fr->Draw();

      canvas->SetEditable(false);

      text_was_set = false;
    }

    void reset()
    {
      for (int cellnum=0; cellnum<nb_cell; ++cellnum)
	{
	  content[cellnum] = 0;
	  content_text_v[cellnum].Clear();
	  cellbox[cellnum].SetFillColor(0);
	}
    }

    float getcontent (int cellnum)
    {
      return content[cellnum];
    }

    void setcontent (int cellnum, float value)
    {
      if (cellnum < nb_cell) content[cellnum] = value;
      else printf("*** wrong cell ID\n");
    }

    void setcontent (int cell_side, int cell_row, int cell_layer, float value)
    {
      int cellnum = cell_side*9*113 + cell_row*9 + cell_layer;
      setcontent(cellnum, value);
    }

    void setcolor (int cell_num, Color_t color)
    {
      if ((cell_num>=0) && (cell_num < 2034))
	{
	  cellbox[cell_num].SetFillColor(color);
	  color_was_set = true;
	}
      else printf("*** wrong cell ID\n");
    }

    void setcolor (int cell_side, int cell_row, int cell_layer, Color_t color)
    {
      int cell_num = cell_side*9*113 + cell_row*9 + cell_layer;
      setcolor(cell_num, color);
    }

    void settext (int cell_num, const char *text)
    {
      if ((cell_num >= 0) && (cell_num < nb_cell))
	{
	  content_text_v[cell_num].SetText(content_text_v[cell_num].GetX(), content_text_v[cell_num].GetY(), text);
	  text_was_set = true;
	}
      else printf("*** wrong cell ID\n");
    }

    void settext (int cell_side, int cell_row, int cell_layer, const char *text)
    {
      int cell_num = cell_side*9*113 + cell_row*9 + cell_layer;
      settext(cell_num, text);
    }

    void fill (int cellnum, float value=1)
    {
      setcontent(cellnum, content[cellnum]+value);
    }

    void update_canvas ()
    {
      canvas->cd();
      canvas->Modified();
      canvas->Update();

      gSystem->ProcessEvents();
    }

    void update (bool update_canvas_too=true)
    {
      if (color_was_set)
	{
	  color_was_set = false;
	  return;
	}

      // autoset Z range [0, content_max]
      // unless setrange() has been called

      float content_min = content[0];
      float content_max = content[0];

      for (int cellnum=1; cellnum<nb_cell; ++cellnum)
	{
	  if (content[cellnum] < content_min) content_min = content[cellnum];
	  if (content[cellnum] > content_max) content_max = content[cellnum];
	}

      content_min = 0;
      if (range_min != -1) content_min = range_min;
      if (range_max != -1) content_max = range_max;

      for (int cellnum=0; cellnum<nb_cell; ++cellnum)
	{
	  if (content[cellnum] != 0)
	    {
	      int color_index = floor (99*(content[cellnum]-content_min)/(content_max-content_min));
	      if (color_index < 0) color_index = 0;
	      else if (color_index >= 100) color_index = 99;
	      cellbox[cellnum].SetFillColor(palette::get_index() + color_index);
	    }
	  else
	    cellbox[cellnum].SetFillColor(0);
	}

      if (palette_axis)
	{
	  palette_histo->SetMinimum(content_min);
	  palette_histo->SetMaximum(content_max);
	  palette_histo->SetContour(100);
	}

      if (update_canvas_too)
	update_canvas();
    }

    TString tracker_name;

    // draw options
    bool draw_cellid;
    bool draw_cellnum;
    bool draw_content;
    TString draw_content_format;
    float range_min, range_max;

    bool color_was_set;
    bool text_was_set;

    TCanvas *canvas;
    TCanvas *canvas_fr;

    TText *label_it;
    TText *label_fr;
    std::vector<TText> row_text_v;

    TH2D *palette_histo;
    TPaletteAxis *palette_axis;

    std::vector<float> content;

    std::vector<TBox>  cellbox;
    std::vector<TText> cellid_text_v;
    std::vector<TText> cellnum_text_v;
    std::vector<TText> content_text_v;

  }; // sndisplay::tracker class

  /////////////////////////////
  // sndisplay::demonstrator //
  /////////////////////////////

  class demonstrator
  {
  public:
    demonstrator (const char *n = "") : demonstrator_name (n)
    {
      canvas = nullptr;

      range_min = range_max = -1;

      // TOP_VIEW //

      const double spacerx = 0.005;
      const double spacery = 0.025;

      const double title_sizey = 0.0615;

      const double mw_sizey = (1-2*spacery-title_sizey)/(2.0 + 4*1.035 + 0.125);
      const double xw_sizey = 1.035*mw_sizey;
      const double se_sizey = 0.125*mw_sizey;
      const double gg_sizey = (1-2*spacery-title_sizey-2*mw_sizey-se_sizey)/18.0;

      const double mw_sizex = (1-2*spacerx)/(20 + 2*0.5*0.720);
      const double xw_sizex = (1-2*spacerx-20*mw_sizex);
      const double se_sizex = (1-2*spacerx-2*xw_sizex);
      const double gg_sizex = se_sizex/113.0;

      // printf("gg_sizex = %f\n", gg_sizex);
      // printf("gg_sizey = %f\n", gg_sizey);

      // MW (column only)

      for (int mw_side=0; mw_side<2; ++mw_side) {

	for (int mw_column=0; mw_column<20; ++mw_column) {

	  double x1 = spacerx + 0.5*xw_sizex + mw_column*mw_sizex;
	  double y1 = spacery + (1-mw_side)*(mw_sizey+4*xw_sizey+se_sizey);

	  double x2 = x1 + mw_sizex;

	  double y2 = y1 + mw_sizey;

	  top_om_content.push_back(0);

	  TBox *box = new TBox(x1, y1, x2, y2);
	  box->SetFillColor(0);
	  box->SetLineWidth(1);
	  top_om_box.push_back(box);

	  TString omid_string = Form("M:%1d.%d.*", mw_side, mw_column);
	  TText *omid_text = new TText (x1+0.5*mw_sizex, y1+0.667*mw_sizey, omid_string);
	  omid_text->SetTextSize(0.032);
	  omid_text->SetTextAlign(22);
	  top_om_text.push_back(omid_text);

	  // TText *content_text = new TText (x1+0.5*mw_sizex, y1+0.333*mw_sizey, "");
	  // content_text->SetTextSize(0.02);
	  // content_text->SetTextAlign(22);
	  // content_text_v.push_back(content_text);

	} // for mw_column

      } // for mw_side

      // XW (column only)

      for (int xw_side=0; xw_side<2; ++xw_side) {

	for (int xw_wall=0; xw_wall<2; ++xw_wall) {

	  for (int xw_column=0; xw_column<2; ++xw_column) {

	    double x1 = spacerx + xw_wall*(xw_sizex+113*gg_sizex);
	    double x2 = x1 + xw_sizex;

	    double y1 = spacery + mw_sizey;

	    if (xw_side == 0)
	      y1 += 2*xw_sizey + se_sizey + xw_column*xw_sizey;
	    else y1 += (1-xw_column)*xw_sizey;

	    double y2 = y1 + xw_sizey;

	    top_om_content.push_back(0);

	    TBox *box = new TBox(x1, y1, x2, y2);
	    box->SetFillColor(0);
	    box->SetLineWidth(1);
	    top_om_box.push_back(box);

	    TString omid_string = Form("X:%1d.%1d.%1d.*", xw_side, xw_wall, xw_column);
	    TText *omid_text = new TText (x1+0.5*xw_sizex, y1+0.6*xw_sizey, omid_string);
	    omid_text->SetTextSize(0.032);
	    omid_text->SetTextAlign(22);
	    top_om_text.push_back(omid_text);

	  }
	}
      }

      for (int gg_side=0; gg_side<2; ++gg_side) {

	  for (int gg_row=0; gg_row<113; ++gg_row) {

	    for (int gg_layer=0; gg_layer<9; ++gg_layer)
	      {
		double x1 = spacerx + xw_sizex + gg_row*gg_sizex;
		double y1 = spacery + mw_sizey;

		if (gg_side == 0)
		  y1 += 9*gg_sizey + se_sizey + gg_layer*gg_sizey;
		else
		  y1 += (8-gg_layer)*gg_sizey;

		double x2 = x1 + gg_sizex;
		double y2 = y1 + gg_sizey;

		top_gg_content.push_back(0);

		TBox *box = new TBox(x1, y1, x2, y2);
		box->SetFillColor(0);
		box->SetLineWidth(1);
		top_gg_box.push_back(box);

		TEllipse *ellipse = new TEllipse((x1+x2)/2, (y1+y2)/2, gg_sizex/2, gg_sizey/2);
		ellipse->SetFillColor(0);
		ellipse->SetLineWidth(1);
		top_gg_ellipse.push_back(ellipse);
	      }
	  }
      }

      title = new TText (spacerx, 1-title_sizey*3/4, "");
      title->SetTextSize(0.056);
      title->SetTextAlign(12);

    } // demonstrator ()


    void setrange(float zmin, float zmax)
    {
      range_min = zmin; range_max = zmax;
    }


    void draw_top()
    {
      update(false);

      if (canvas == nullptr)
	{
	  const int canvas_width  = 1600;
	  const int canvas_height = 400; // 374 without title

	  canvas = new TCanvas (Form("C_demonstrator_%s",demonstrator_name.Data()), Form("%s",demonstrator_name.Data()), canvas_width, canvas_height);

	  // force canvas exact size
	  int decoration_width = canvas_width - canvas->GetWw();
	  int decoration_height = canvas_height - canvas->GetWh();
	  canvas->SetWindowSize(canvas_width+decoration_width, canvas_height+decoration_height);

	  // preserve width/height ratio in case of resizing
	  canvas->SetFixedAspectRatio();
	}
      else canvas->cd();

      for (int mw_side=0; mw_side<2; ++mw_side)
	{
	  for (int mw_column=0; mw_column<20; ++mw_column)
	    {
	      int top_om_num = mw_side*20 + mw_column;
 	      top_om_box[top_om_num]->Draw("l");
	      top_om_text[top_om_num]->Draw();
	    }
	}

      for (int xw_side=0; xw_side<2; ++xw_side)
	{
	  for (int xw_wall=0; xw_wall<2; ++xw_wall)
	    {
	      for (int xw_column=0; xw_column<2; ++xw_column)
		{
		  int top_om_num = 40 + xw_side*2*2 + xw_wall*2 + xw_column;
		  top_om_box[top_om_num]->Draw("l");
		  top_om_text[top_om_num]->Draw();
		}
	    }
	}

      for (int gg_side=0; gg_side<2; ++gg_side)
	{
	  for (int gg_row=0; gg_row<113; ++gg_row)
	    {
	      for (int gg_layer=0; gg_layer<9; ++gg_layer)
		{
		  int top_gg_num = gg_side*113*9 + gg_row*9 + gg_layer;
		  top_gg_box[top_gg_num]->Draw("l");
		  top_gg_ellipse[top_gg_num]->Draw("l");
		  // top_g_text[top_gg_num]->Draw();
		}
	    }
	}

      title->Draw();

    } // draw_top

    void setomcontent (int om_num, float value)
    {
      int top_om_num = -1;

      if (om_num < 260) // MW IT
	{
	  int om_side = 0;
	  int om_column = (om_num/13);
	  top_om_num = om_side*20 + om_column;
	}
      else if (om_num < 520) // MW IT
	{
	  int om_side = 1;
	  int om_column = (om_num-260)/13;
	  top_om_num = om_side*20 + om_column;
	}
      else if (om_num < 648) // XW
	{
	  int om_side = (om_num < 584) ? 0 : 1;
	  int om_wall = (om_num-520-om_side*64)/32;
	  int om_column = (om_num-520-om_side*64-om_wall*32)/16;
	  top_om_num = 40 + om_side*2*2 + om_wall*2 + om_column;
	}

      top_om_content[top_om_num] = value;
    }


    void setggcontent (int cell_num, float value)
    {
      if (cell_num < 2034) top_gg_content[cell_num] = value;
      else printf("*** wrong cell ID\n");
    }

    float getggcontent (int cell_num)
    {
      return top_gg_content[cell_num];
    }

    void setggcontent (int cell_side, int cell_row, int cell_layer, float value)
    {
      int cell_num = cell_side*9*113 + cell_row*9 + cell_layer;
      setggcontent(cell_num, value);
    }

    void setggcolor (int cell_num, Color_t color)
    {
      if (cell_num < 2034) top_gg_ellipse[cell_num]->SetFillColor(color);
      else printf("*** wrong cell ID\n");
    }

    void setggcolor (int cell_side, int cell_row, int cell_layer, Color_t color)
    {
      int cell_num = cell_side*9*113 + cell_row*9 + cell_layer;
      setggcolor(cell_num, color);
    }

    void settitle (const char *text)
    {
      title->SetText(title->GetX(), title->GetY(), text);
    }

    void reset ()
    {
      for (size_t om=0; om<top_om_content.size(); ++om)
	{
	  top_om_content[om] = 0;
	  top_om_box[om]->SetFillColor(0);
	}

      for (size_t gg=0; gg<top_gg_content.size(); ++gg)
	{
	  top_gg_content[gg] = 0;
	  top_gg_ellipse[gg]->SetFillColor(0);
	  // top_gg_box[gg]->SetFillColor(0);
	}
    }

    void update_canvas ()
    {
      canvas->Modified();
      canvas->Update();

      gSystem->ProcessEvents();
    }

    void update (bool update_canvas_too=true)
    {
      float top_content_min = top_om_content[0];
      float top_content_max = top_om_content[0];

      for (size_t om=0; om<top_om_content.size(); ++om)
    	{
    	  if (top_om_content[om] < top_content_min) top_content_min = top_om_content[om];
    	  if (top_om_content[om] > top_content_max) top_content_max = top_om_content[om];
    	}

      for (size_t gg=0; gg<top_gg_content.size(); ++gg)
    	{
    	  if (top_gg_content[gg] < top_content_min) top_content_min = top_gg_content[gg];
    	  if (top_gg_content[gg] > top_content_max) top_content_max = top_gg_content[gg];
    	}

      top_content_min = 0;
      if (range_min != -1) top_content_min = range_min;
      if (range_max != -1) top_content_max = range_max;
      // printf("Z range = [%f, %f] for '%s'\n", top_content_min, top_content_max, demonstrator_name.Data());

      for (size_t om=0; om<top_om_content.size(); ++om)
    	{
    	  if (top_om_content[om] != 0)
    	    {
    	      int color_index = floor (99*(top_om_content[om]-top_content_min)/(top_content_max-top_content_min));
    	      if (color_index < 0) color_index = 0;
    	      else if (color_index >= 100) color_index = 99;
	      top_om_box[om]->SetFillColor(palette::get_index() + color_index);
    	    }
    	  else
    	    top_om_box[om]->SetFillColor(0);
    	}

      for (size_t gg=0; gg<top_gg_content.size(); ++gg)
    	{
    	  if (top_gg_content[gg] != 0)
    	    {
    	      int color_index = floor (99*(top_gg_content[gg]-top_content_min)/(top_content_max-top_content_min));
    	      if (color_index < 0) color_index = 0;
    	      else if (color_index >= 100) color_index = 99;
	      top_gg_ellipse[gg]->SetFillColor(palette::get_index() + color_index);
    	    }
    	  else
    	    top_gg_ellipse[gg]->SetFillColor(0);
    	}

      if (update_canvas_too)
	update_canvas();
    }

    //

    TString demonstrator_name;

    TCanvas *canvas;
    TText *title;

    float range_min, range_max;

    std::vector<float> top_om_content;
    std::vector<TBox*> top_om_box;
    std::vector<TText*> top_om_text;

    std::vector<float> top_gg_content;
    std::vector<TBox*> top_gg_box;
    std::vector<TEllipse*> top_gg_ellipse;
    // std::vector<TText*>  top_gg_text;

  }; // sndisplay::demonstrator class

} // sndisplay namespace

////////////////////////////////////////////////////////////////

// sndisplay_calorimeter_test_values()
// => sndisplay::calorimeter usage by filling OMs content with random value

sndisplay::calorimeter *sncalo = nullptr;

void sndisplay_calorimeter_test_values (bool with_palette = true)
{
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);

  sncalo->draw_content_label("%.0f");

  int ncalo_tot;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;

  TFile *file = new TFile("data/simu_Bi_207.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_number);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);


  double rate[712];
  for (int i = 0; i < tree->GetEntries()/10; i++) {
    tree->GetEntry(i);
    // cout << "size = " << om_number->size() << endl;
    for (int j = 0; j < om_number->size(); j++) {
      // if (energyvis_ubc->at(j) > 0.6) {
        rate[om_number->at(j)]++;
        // cout << "ok" << endl;
      // }
    }
  }


  for (int omnum=0; omnum<520; ++omnum) // MW
    sncalo->setcontent(omnum, rate[omnum]);

  for (int omnum=520; omnum<648; ++omnum) // XW
    sncalo->setcontent(omnum, 0);

  for (int omnum=648; omnum<712; ++omnum) // GV
    sncalo->setcontent(omnum, 0);

  sncalo->setrange(0, 50000); // to force z axis range

  sncalo->draw();

  sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}
//
void sndisplay_calorimeter_cut_rate_e ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);

  sncalo->draw_content_label("%.1f");

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *charge = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<int> *flag = new std::vector<int>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;

  TFile *file = new TFile("cut_bi_10.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &om_number);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress("flag", &flag);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);

  TH2D* Gamma = new TH2D("Gamma","Gamma",520,0,520,1000,0,2);

  int e_compteur;
  int compteur =0;

  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    e_compteur = 0;
    // for (int j = 0; j < flag_e_event->size(); j++)
    //   if (flag_e_event->at(j) == 1 && energyvis_ubc->at(j)>0.7 && ncalo_tot > 1) e_compteur++;

    for (int j = 0; j < flag_e_event->size(); j++) {
      if (flag_e_event->at(j) == 1 && energyvis_ubc->at(j) > 0.7 && ncalo_tot <2) {
      // if (flag_e_event->at(j) == 0 && e_compteur > 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {

        rate_simu[om_number->at(j)]++;
        compteur++;
      }
    }
  }

  cout << "total simu event : " << compteur << endl;
  std::vector<long> *timestamp = new std::vector<long>;

  TFile *file2 = new TFile("../cut_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("flag",1);
  tree2->SetBranchAddress("flag", &flag);
  tree2->SetBranchStatus("flag_e_event",1);
  tree2->SetBranchAddress("flag_e_event", &flag_e_event);
  tree2->SetBranchStatus("flag_charge",1);
  tree2->SetBranchAddress("flag_charge", &flag_charge);
  tree2->SetBranchStatus("flag_associated_nohit",1);
  tree2->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree2->SetBranchStatus("flag_MW",1);
  tree2->SetBranchAddress("flag_MW", &flag_MW);
  tree2->SetBranchStatus("flag_calo_square",1);
  tree2->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree2->SetBranchStatus("flag_source_square",1);
  tree2->SetBranchAddress("flag_source_square", &flag_source_square);
  tree2->SetBranchStatus("flag_last_column",1);
  tree2->SetBranchAddress("flag_last_column", &flag_last_column);
  tree2->SetBranchStatus("flag_last_z",1);
  tree2->SetBranchAddress("flag_last_z", &flag_last_z);
  tree2->SetBranchStatus("charge",1);
  tree2->SetBranchAddress("charge", &charge);
  tree2->SetBranchStatus("z_last_gg",1);
  tree2->SetBranchAddress("z_last_gg", &z_last_gg);
  tree2->SetBranchStatus("calo_nohit_om_time",1);
  tree2->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);
  tree2->SetBranchStatus("energy",1);
  tree2->SetBranchAddress("energy", &energyvis_ubc);
  tree2->SetBranchStatus("calo_timestamp",1);
  tree2->SetBranchAddress("calo_timestamp", &timestamp);

  int last_j;
  double rate[712];
  memset (rate, 0, 712*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    e_compteur = 0;

    // for (int j = 0; j < flag_e_event->size(); j++)
    //   if (flag_e_event->at(j) == 1 ) {
    //     e_compteur++;
    //     last_j = j;
    //   }
    // if (e_compteur <1) continue;
    for (int j = 0; j < om_number->size(); j++) {

      // float time = abs(timestamp->at(j)-timestamp->at(last_j));
            // cout << time << endl;/
      // if (time < 16) cout << time << " and "
      // if ( energyvis_ubc->at(j)>0.7 && calo_nohit_om_time->at(j) < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1 && flag_last_z->at(j) == 1)  {
      // if (time < 16 && flag_e_event->at(j) == 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {
      if (energyvis_ubc->at(j) > 0.7 && flag_e_event->at(j) == 1 && calo_nohit_om_time->at(j) < 2) {
              // if (om_number->at(j) == 175 && energyvis_ubc->at(j) < 0.4 && energyvis_ubc->at(j) > 0.38) cout << "entry : " << i << endl;
        // Gamma->Fill(om_number->at(j), energyvis_ubc->at(j));
        rate[om_number->at(j)]++;
        compteur2++;
      }
      // time = 0;
    }
  }
  // TFile *newfile = new TFile("test.root", "RECREATE");
  // newfile->cd();
  // Gamma->Write();
  // newfile->Close();
  // return;
  //
  cout << "total data event : " << compteur2 << endl;
  // for (size_t i = 0; i < 520; i++) {
  //   std::cout << "rate simu = " << rate_simu[i] << " rate = " << rate[i] << '\n';
  // }



  for (int omnum=0; omnum<520; ++omnum) // MW
    if (rate_simu[omnum] > 0 && rate[omnum] > 0)
    sncalo->setcontent(omnum, rate[omnum]);
    // sncalo->setcontent(omnum, rate_simu[omnum]/(rate[omnum]/1.78));

  sncalo->setrange(-1, -1); // to force z axis range

  sncalo->draw();
  sncalo->draw1();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}

void sndisplay_calorimeter_cut_rate_gamma ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");


  double simu_mean[520];
  memset (simu_mean, 0, 520*sizeof(double));
  double correction_tab[520];
  memset (correction_tab, 0, 520*sizeof(double));
  int om_num;
  double mean, correction;

  TFile *correction_file = new TFile("Bi_fit/ubc_corrector.root", "READ");
  TTree* correction_tree = (TTree*)correction_file->Get("Result_tree");
  correction_tree->SetBranchStatus("*",0);
  correction_tree->SetBranchStatus("om_number",1);
  correction_tree->SetBranchAddress("om_number", &om_num);
  correction_tree->SetBranchStatus("mean",1);
  correction_tree->SetBranchAddress("mean", &mean);
  correction_tree->SetBranchStatus("correction",1);
  correction_tree->SetBranchAddress("correction", &correction);

  for (int i = 0; i < correction_tree->GetEntries(); i++) {
    correction_tree->GetEntry(i);
    simu_mean[om_num] = mean;
    correction_tab[om_num] = correction;
  }
  cout << "ubc corrector ok" << endl;

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis = new std::vector<double>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *charge = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<int> *flag = new std::vector<int>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;

  TFile *file = new TFile("cut_bi_new.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &om_number);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress("flag", &flag);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
  tree->SetBranchStatus("energyvis",1);
  tree->SetBranchAddress("energyvis", &energyvis);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);

  TH2D* Gamma_data = new TH2D("Gamma_data","Gamma_data",520,0,520,100,0,2);
  TH2D* Gamma_simu_ubc = new TH2D("Gamma_simu_ubc","Gamma_simu_ubc",520,0,520,100,0,2);
  TH2D* Gamma_simu = new TH2D("Gamma_simu","Gamma_simu",520,0,520,100,0,2);
  TH1D* simu = new TH1D("simu","simu",10,0,10);
  int e_compteur;
  int compteur =0;

  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries()/10; i++) {
    tree->GetEntry(i);
    e_compteur = 0;
    for (int j = 0; j < flag_e_event->size(); j++)
      if (flag_e_event->at(j) == 1 && energyvis_ubc->at(j)>0.7 && ncalo_tot > 1) {
        Gamma_simu_ubc->Fill(om_number->at(j), energyvis_ubc->at(j)*correction_tab[om_number->at(j)]);
        Gamma_simu->Fill(om_number->at(j), energyvis->at(j));

        e_compteur++;
      }
    if (e_compteur <1) continue;
    for (int j = 0; j < flag_e_event->size(); j++) {
      if (flag_e_event->at(j) == 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {
        // if (flag_e_event->at(j) == 1 && energyvis_ubc->at(j) > 0.7 && ncalo_tot <2) {
        // Gamma_simu_ubc->Fill(om_number->at(j), energyvis_ubc->at(j)*correction_tab[om_number->at(j)]);
        // Gamma_simu->Fill(om_number->at(j), energyvis->at(j));
        simu->Fill(ncalo_tot);
        rate_simu[om_number->at(j)]++;
        compteur++;
      }
    }
  }

  cout << "total simu event : " << compteur << endl;
  std::vector<long> *timestamp = new std::vector<long>;

  TFile *file2 = new TFile("../cut_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("flag",1);
  tree2->SetBranchAddress("flag", &flag);
  tree2->SetBranchStatus("flag_e_event",1);
  tree2->SetBranchAddress("flag_e_event", &flag_e_event);
  tree2->SetBranchStatus("flag_charge",1);
  tree2->SetBranchAddress("flag_charge", &flag_charge);
  tree2->SetBranchStatus("flag_associated_nohit",1);
  tree2->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree2->SetBranchStatus("flag_MW",1);
  tree2->SetBranchAddress("flag_MW", &flag_MW);
  tree2->SetBranchStatus("flag_calo_square",1);
  tree2->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree2->SetBranchStatus("flag_source_square",1);
  tree2->SetBranchAddress("flag_source_square", &flag_source_square);
  tree2->SetBranchStatus("flag_last_column",1);
  tree2->SetBranchAddress("flag_last_column", &flag_last_column);
  tree2->SetBranchStatus("flag_last_z",1);
  tree2->SetBranchAddress("flag_last_z", &flag_last_z);
  tree2->SetBranchStatus("charge",1);
  tree2->SetBranchAddress("charge", &charge);
  tree2->SetBranchStatus("z_last_gg",1);
  tree2->SetBranchAddress("z_last_gg", &z_last_gg);
  tree2->SetBranchStatus("calo_nohit_om_time",1);
  tree2->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);
  tree2->SetBranchStatus("energy",1);
  tree2->SetBranchAddress("energy", &energyvis_ubc);
  tree2->SetBranchStatus("calo_timestamp",1);
  tree2->SetBranchAddress("calo_timestamp", &timestamp);

  TH1D* data = new TH1D("data","data",10,0,10);
  int last_j;
  double rate[712];
  memset (rate, 0, 712*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries()/10; i++) {
    tree2->GetEntry(i);
    e_compteur = 0;

    for (int j = 0; j < flag_e_event->size(); j++)
      if (flag_e_event->at(j) == 1  && energyvis_ubc->at(j)>0.7 && calo_nohit_om_time->at(j)>1) {
        Gamma_data->Fill(om_number->at(j), energyvis_ubc->at(j)*simu_mean[om_number->at(j)]);
        // if (om_number->at(j) == 175) cout << "entry : " << i << " -> om:" << om_number->at(j) << " energy :" << energyvis_ubc->at(j) << endl;
        e_compteur++;
        last_j = j;
      }
    if (e_compteur <1) continue;

    for (int j = 0; j < om_number->size(); j++) {

      float time = abs(timestamp->at(j)-timestamp->at(last_j));
      if (time < 16 && flag_e_event->at(j) == 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {
        rate[om_number->at(j)]++;
            data->Fill(calo_nohit_om_time->at(last_j)+1);
        compteur2++;
      }
      time = 0;
    }
  }
  simu->Draw();
  data->Draw("same");
  simu->SetLineColor(kRed);
  return;
  TFile *newfile = new TFile("test_e.root", "RECREATE");
  newfile->cd();
  Gamma_data->Write();
  Gamma_simu->Write();
  Gamma_simu_ubc->Write();
  newfile->Close();

  //
  cout << "total data event : " << compteur2 << endl;




  for (int omnum=0; omnum<520; ++omnum) // MW
    if (rate_simu[omnum] > 0 && rate[omnum] > 0)
    sncalo->setcontent(omnum, rate_simu[omnum]/(rate[omnum]/1.78));
    // sncalo->setcontent(omnum, rate[omnum]);

  sncalo->setrange(-1, -1); // to force z axis range

  sncalo->draw();
  sncalo->draw1();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

}

void sndisplay_calorimeter_sources ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);

  sncalo->draw_content_label("%.1f");

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *charge = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<double> *z_last_gg = new std::vector<double>;
  std::vector<int> *flag = new std::vector<int>;
  std::vector<int> *calo_nohit_om_time = new std::vector<int>;
  std::vector<int> *source_number = new std::vector<int>;

  TFile *file = new TFile("cut_bi_good_gas.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_num",1);
  tree->SetBranchAddress("om_num", &om_number);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress("flag", &flag);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energyvis_ubc",1);
  tree->SetBranchAddress("energyvis_ubc", &energyvis_ubc);
  tree->SetBranchStatus("z_last_gg",1);
  tree->SetBranchAddress("z_last_gg", &z_last_gg);
  tree->SetBranchStatus("source_number",1);
  tree->SetBranchAddress("source_number", &source_number);

  TH2D* Gamma = new TH2D("Gamma","Gamma",520,0,520,1000,0,2);

  int e_compteur;
  int compteur =0;

  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    for (int j = 0; j < flag_e_event->size(); j++) {
      if (flag_e_event->at(j) == 1 && energyvis_ubc->at(j) > 0.7 ) {
      // if (flag_e_event->at(j) == 0 && e_compteur > 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {

        rate_simu[om_number->at(j)]++;
        compteur++;
      }
    }
  }

  cout << "total simu event : " << compteur << endl;
  // return;
  std::vector<long> *timestamp = new std::vector<long>;

  TFile *file2 = new TFile("../cut_974.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("flag",1);
  tree2->SetBranchAddress("flag", &flag);
  tree2->SetBranchStatus("flag_e_event",1);
  tree2->SetBranchAddress("flag_e_event", &flag_e_event);
  tree2->SetBranchStatus("flag_charge",1);
  tree2->SetBranchAddress("flag_charge", &flag_charge);
  tree2->SetBranchStatus("flag_associated_nohit",1);
  tree2->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree2->SetBranchStatus("flag_MW",1);
  tree2->SetBranchAddress("flag_MW", &flag_MW);
  tree2->SetBranchStatus("flag_calo_square",1);
  tree2->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree2->SetBranchStatus("flag_source_square",1);
  tree2->SetBranchAddress("flag_source_square", &flag_source_square);
  tree2->SetBranchStatus("flag_last_column",1);
  tree2->SetBranchAddress("flag_last_column", &flag_last_column);
  tree2->SetBranchStatus("flag_last_z",1);
  tree2->SetBranchAddress("flag_last_z", &flag_last_z);
  tree2->SetBranchStatus("charge",1);
  tree2->SetBranchAddress("charge", &charge);
  tree2->SetBranchStatus("z_last_gg",1);
  tree2->SetBranchAddress("z_last_gg", &z_last_gg);
  tree2->SetBranchStatus("calo_nohit_om_time",1);
  tree2->SetBranchAddress("calo_nohit_om_time", &calo_nohit_om_time);
  tree2->SetBranchStatus("energy",1);
  tree2->SetBranchAddress("energy", &energyvis_ubc);
  tree2->SetBranchStatus("calo_timestamp",1);
  tree2->SetBranchAddress("calo_timestamp", &timestamp);
  tree2->SetBranchStatus("source_number",1);
  tree2->SetBranchAddress("source_number", &source_number);

  int last_j;
  double rate[712];
  memset (rate, 0, 712*sizeof(double));
  int compteur2 =0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    e_compteur = 0;

    for (int j = 0; j < om_number->size(); j++) {

      // float time = abs(timestamp->at(j)-timestamp->at(last_j));

      // if (time < 16) cout << time << " and "
      if ( energyvis_ubc->at(j)>0.7  && flag_e_event->at(j) == 1)  {
      // if (time < 16 && flag_e_event->at(j) == 0 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_calo_square->at(j) == 0 && energyvis_ubc->at(j) > 0.22) {
      // if (energyvis_ubc->at(j) > 0.7 && flag_e_event->at(j) == 1 && ) {
        // cout << "i: " << i << " : source_number " << source_number->at(j) << " : energy " << energyvis_ubc->at(j)<< endl;
        rate[om_number->at(j)]++;
        compteur2++;
      }
      // time = 0;
    }
  }

  cout << "total data event : " << compteur2 << endl;

  int om;
  float teres;

  TFile *file3 = new TFile("Eres.root", "READ");
  TTree* tree3 = (TTree*)file3->Get("Result_tree");
  tree3->SetBranchStatus("*",0);
  tree3->SetBranchStatus("eres",1);
  tree3->SetBranchAddress("eres", &teres);
  tree3->SetBranchStatus("om",1);
  tree3->SetBranchAddress("om", &om);
  std::cout << "/* message */" << '\n';

  double eres[520];
  memset (eres, 0, 520*sizeof(double));
  for (int i = 0; i < tree3->GetEntries(); i++) {
    tree3->GetEntry(i);
    eres[om] = teres;
  }

  TH2F* rate_eres = new TH2F("rate_eres","rate_eres", 200, 0, 2, 150,10,25);
  TH2F* rate_eres_bad = new TH2F("rate_eres_bad_it_zone","rate_eres_bad_it_zone", 200, 0, 2, 150,10,25);
  TH2F* rate_eres_bad2 = new TH2F("rate_eres_bad_it_other","rate_eres_it_bad_other", 200, 0, 2, 150,10,25);
  TH2F* rate_eres_bad3 = new TH2F("rate_eres_bad_fr_zone","rate_eres_bad_fr_zone", 200, 0, 2, 150,10,25);
  TH2F* rate_eres_bad4 = new TH2F("rate_eres_bad_fr_other","rate_eres_bad_fr_other", 200, 0, 2, 150,10,25);

  int bad[9] = {118,120,121,122,123,131,132,133,134};
  int bad2[4] = {2,5,84,192};
  int bad3[8] = {508,510,511,512, 513,514,515,516};
  int bad4[13] = {292,316,330,340,348,366, 368,370,458,470,471,491,495};
  int c = 0;
  int c2 = 0;
  int c3 = 0;
  int c4 = 0;
  for (int omnum=0; omnum<520; ++omnum){ // MW
    if (rate_simu[omnum] > 0 && rate[omnum] > 0){
      sncalo->setcontent(omnum, rate_simu[omnum]/(rate[omnum]/1.78));
      if (omnum == bad[c]) {
        rate_eres_bad->Fill(rate_simu[omnum]/(rate[omnum]/1.78), eres[omnum]);
        c++;
      }
      else if (omnum == bad2[c2]) {
        rate_eres_bad2->Fill(rate_simu[omnum]/(rate[omnum]/1.78), eres[omnum]);
        c2++;
      }
      else if (omnum == bad3[c3]) {
        rate_eres_bad3->Fill(rate_simu[omnum]/(rate[omnum]/1.78), eres[omnum]);
        c3++;
      }
      else if (omnum == bad4[c4]) {
        rate_eres_bad4->Fill(rate_simu[omnum]/(rate[omnum]/1.78), eres[omnum]);
        c4++;
      }
      else if (omnum%13 !=1 && omnum%13!=11) {
        rate_eres->Fill(rate_simu[omnum]/(rate[omnum]/1.78), eres[omnum]);
      }

    }
  }
  rate_eres->GetXaxis()->SetTitle("Rate simu/data");
  rate_eres->GetYaxis()->SetTitle("Eres");
    // sncalo->setcontent(omnum, rate[omnum]);
  TLegend *legend = new TLegend(0.5,0.7,0.7,0.9);
  legend->AddEntry(rate_eres, "rate/eres normal OM", "P");
  legend->AddEntry(rate_eres_bad, "rate/eres bad OM it zone", "P");
  legend->AddEntry(rate_eres_bad2, "rate/eres bad OM it other", "P");
  legend->AddEntry(rate_eres_bad3, "rate/eres bad OM fr zone", "P");
  legend->AddEntry(rate_eres_bad4, "rate/eres bad OM fr other", "P");

  rate_eres->Draw();
  rate_eres_bad->Draw("same");
  rate_eres_bad2->Draw("same");
  rate_eres_bad3->Draw("same");
  rate_eres_bad4->Draw("same");
  legend->Draw("same");

  rate_eres->SetMarkerStyle(2);
  rate_eres_bad->SetMarkerColor(kRed);
  rate_eres_bad->SetMarkerStyle(8);
  rate_eres_bad2->SetMarkerColor(kRed);
  rate_eres_bad2->SetMarkerStyle(4);
  rate_eres_bad3->SetMarkerColor(kBlue);
  rate_eres_bad3->SetMarkerStyle(22);
  rate_eres_bad4->SetMarkerColor(kBlue);
  rate_eres_bad4->SetMarkerStyle(26);
  sncalo->setrange(-1, -1); // to force z axis range

  sncalo->draw();
  sncalo->draw1();

}

void sndisplay_calorimeter_simu_cut_rate ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);

  sncalo->draw_content_label("%.3f");

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;

  TFile *file = new TFile("../cut_974_no.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);

  int compteur =0;
  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries()/10; i++) {
    tree->GetEntry(i);
    for (int j = 0; j < flag_e_event->size(); j++) {
      // if (energyvis_ubc->at(j) > 0.6 && ncalo_tot < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1)  {
      if (ncalo_tot < 2 && flag_e_event->at(j) == 1 && energy->at(j) > 0.7) {
        rate_simu[om_number->at(j)]++;
        compteur++;
      }
    }
  }
  cout << "total simu event : " << compteur << endl;

  TFile *file2 = new TFile("cut_bi.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("ncalo_tot",1);
  tree2->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree2->SetBranchStatus("om_num",1);
  tree2->SetBranchAddress("om_num", &om_number);
  tree2->SetBranchStatus("flag_e_event",1);
  tree2->SetBranchAddress("flag_e_event", &flag_e_event);
  tree2->SetBranchStatus("flag_charge",1);
  tree2->SetBranchAddress("flag_charge", &flag_charge);
  tree2->SetBranchStatus("flag_associated_nohit",1);
  tree2->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree2->SetBranchStatus("flag_MW",1);
  tree2->SetBranchAddress("flag_MW", &flag_MW);
  tree2->SetBranchStatus("flag_calo_square",1);
  tree2->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree2->SetBranchStatus("flag_source_square",1);
  tree2->SetBranchAddress("flag_source_square", &flag_source_square);
  tree2->SetBranchStatus("flag_last_column",1);
  tree2->SetBranchAddress("flag_last_column", &flag_last_column);
  tree2->SetBranchStatus("flag_last_z",1);
  tree2->SetBranchAddress("flag_last_z", &flag_last_z);
  tree2->SetBranchStatus("energyvis_ubc",1);
  tree2->SetBranchAddress("energyvis_ubc", &energyvis_ubc);

  int compteur2 = 0;
  double rate_simu2[712];
  memset (rate_simu2, 0, 712*sizeof(double));
  for (int i = 0; i < tree2->GetEntries()/10; i++) {
    tree2->GetEntry(i);
    for (int j = 0; j < flag_e_event->size(); j++) {
      // if (energyvis_ubc->at(j) > 0.6 && ncalo_tot < 2 && flag_e_event->at(j) == 1) {
      // if (energyvis_ubc->at(j) > 0.6 && ncalo_tot < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1 ) {
      if (ncalo_tot < 2 && flag_e_event->at(j) == 1 && energyvis_ubc->at(j) > 0.7) {
        rate_simu2[om_number->at(j)]++;
        compteur2++;
      }
    }
  }
  cout << "total simu event : " << compteur2 << endl;

  int om;
  float teres;

  TFile *file3 = new TFile("Eres.root", "READ");
  TTree* tree3 = (TTree*)file3->Get("Result_tree");
  tree3->SetBranchStatus("*",0);
  tree3->SetBranchStatus("eres",1);
  tree3->SetBranchAddress("eres", &teres);
  tree3->SetBranchStatus("om",1);
  tree3->SetBranchAddress("om", &om);

  double eres[520];
  memset (eres, 0, 520*sizeof(double));
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    eres[om] = teres;
  }

  //
  // cout << "total data event : " << compteur2 << endl;
  // for (size_t i = 0; i < 520; i++) {
  //   std::cout <<"om : " << i << "  rate simu dead = " << rate_simu[i] << " rate simu perfect= " << rate_simu2[i] << '\n';
  // }

  TH2F* rate_eres = new TH2F("rate_eres","rate_eres", 200, 0, 2, 150,10,25);

  for (int omnum=0; omnum<520; ++omnum) // MW
    if (rate_simu[omnum] > 0){
     sncalo->setcontent(omnum, rate_simu2[omnum]/(rate_simu[omnum]/1.78));
     // sncalo->setcontent(omnum, rate_simu[omnum]);
     rate_eres->Fill(rate_simu2[omnum]/(rate_simu[omnum]), eres[omnum]);
   }
  // sncalo->setrange(0, 3); // to force z axis range
  rate_eres->Draw();
  sncalo->draw();
  sncalo->draw1();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}

void sndisplay_calorimeter_taux_source ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);

  sncalo->draw_content_label("%.0f");

  int ncalo_tot;
  std::vector<int> *flag_e_event = new std::vector<int>;
  std::vector<int> *om_number = new std::vector<int>;
  std::vector<double> *energyvis_ubc = new std::vector<double>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<int> *flag_charge = new std::vector<int>;
  std::vector<int> *flag_associated_nohit = new std::vector<int>;
  std::vector<int> *flag_MW = new std::vector<int>;
  std::vector<int> *flag_calo_square = new std::vector<int>;
  std::vector<int> *flag_source_square = new std::vector<int>;
  std::vector<int> *flag_last_column = new std::vector<int>;
  std::vector<int> *flag_last_z = new std::vector<int>;
  std::vector<int> *source_number = new std::vector<int>;

  TFile *file = new TFile("../cut_974.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("ncalo_tot",1);
  tree->SetBranchAddress("ncalo_tot", &ncalo_tot);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("flag_e_event",1);
  tree->SetBranchAddress("flag_e_event", &flag_e_event);
  tree->SetBranchStatus("flag_charge",1);
  tree->SetBranchAddress("flag_charge", &flag_charge);
  tree->SetBranchStatus("flag_associated_nohit",1);
  tree->SetBranchAddress("flag_associated_nohit", &flag_associated_nohit);
  tree->SetBranchStatus("flag_MW",1);
  tree->SetBranchAddress("flag_MW", &flag_MW);
  tree->SetBranchStatus("flag_calo_square",1);
  tree->SetBranchAddress("flag_calo_square", &flag_calo_square);
  tree->SetBranchStatus("flag_source_square",1);
  tree->SetBranchAddress("flag_source_square", &flag_source_square);
  tree->SetBranchStatus("flag_last_column",1);
  tree->SetBranchAddress("flag_last_column", &flag_last_column);
  tree->SetBranchStatus("flag_last_z",1);
  tree->SetBranchAddress("flag_last_z", &flag_last_z);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchStatus("source_number",1);
  tree->SetBranchAddress("source_number", &source_number);

  int compteur =0;
  double rate_simu[712];
  memset (rate_simu, 0, 712*sizeof(double));
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    for (int j = 0; j < flag_e_event->size(); j++) {
      // if (energyvis_ubc->at(j) > 0.6 && ncalo_tot < 2 && flag_charge->at(j) == 1 && flag_MW->at(j) == 1 && flag_associated_nohit->at(j)==1 && flag_last_column->at(j) == 1 && flag_calo_square->at(j) == 1 && flag_source_square->at(j) == 1)  {
      // if (flag_e_event->at(j) == 1 && source_number->at(j) !=204 &&  source_number->at(j) !=201 &&source_number->at(j) !=202 &&source_number->at(j) !=203 &&source_number->at(j)) {
      // if (flag_e_event->at(j) == 1 && (source_number->at(j) ==304 ||  source_number->at(j) ==301 ||source_number->at(j) ==302 ||source_number->at(j) ==303 )) {
      if (flag_e_event->at(j) == 1 && (source_number->at(j) ==304)) {

        rate_simu[om_number->at(j)]++;
        compteur++;
      }
    }
  }
  cout << "total simu event : " << compteur << endl;

  //
  // cout << "total data event : " << compteur2 << endl;
  // for (size_t i = 0; i < 520; i++) {
  //   std::cout <<"om : " << i << "  rate simu dead = " << rate_simu[i] << " rate simu perfect= " << rate_simu2[i] << '\n';
  // }


  for (int omnum=0; omnum<520; ++omnum) // MW
    if (rate_simu[omnum] > 0){
     sncalo->setcontent(omnum, (rate_simu[omnum]));
     // sncalo->setcontent(omnum, rate_simu[omnum]);
   }
  // sncalo->setrange(0, 3); // to force z axis range
  sncalo->draw();
  sncalo->draw1();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}

void sndisplay_calorimeter_taux_rapport()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);

  sncalo->draw_content_label("%.1f");

  double ncal[520] = {0,1/2.7,1/0.8,1/0.9,1/0.7,1/0.8,1/0.8,1/0.9,1/0.8,1/0.9,1/0.9,1/2.6,0,0,1/2.2,1/0.8,1/0.9,1/0.9,1/1.1,1/0.9,1/0.9,1,1/0.9,1/0.9,1/4.,0,0,1/2.2,1/0.8,1/0.8,1/0.9,1/0.8,1/0.9,1/0.9,1,1/0.9,1/0.9,1/3.3,0,0,1/3.0,1/0.9,1/0.8,1/0.8,1/0.8,1/0.9,1/0.8,1/0.9,1,1/0.9,1/4.1,0,0,1/2.3,1/0.8,1/0.8,1/0.9,1/0.8,1/0.8,1/0.9,1/0.8,1/0.8,1/0.9,1/3.0,0,0,1/2.2,1/0.9,1/0.9,1/0.9,1/0.9,1/0.8,1/0.7,1/0.9,0,1/0.9,1/2.6,0,0,1/2.3,0,1/0.7,1/0.7,1/0.7,1/0.8,1/0.8,1/0.9,1/0.8,1/0.8,1/4.9,0,0,1/2.5,1/0.8,1/0.8,1/0.8,1/0.7,1/0.8,0,1/0.9,1/0.9,1/0.9,1/3.7,0,0,1/1.9,1/0.8,1/0.8,1/0.9,1/0.9,1/0.8,1/0.8,1/0.9,1/0.9,1/0.9,1/2.7,0,0,1/4.3,0,1/1.2,1/1.1,1,1/1.1,1/0.9,1/0.9,1/0.9,1,1/4.5,0,0,1/4.,1/1.5,1/1.2,1,1/1.1,0,1,1/0.9,1/0.9,1,1/2.6,0,0,1/1.8,1/0.8,0,1/0.8,1/0.8,1/0.9,1,1/0.8,1/0.9,1/0.9,1/2.6,0,0,1/2.9,1/0.8,1/0.8,1/0.8,1/0.8,1/0.9,1/0.9,1,1/0.9,1/0.9,1/3.5,0,0,1/3.1,1/0.9,1/0.8,1/0.9,1/1.1,1/0.8,1,1/0.8,1/0.7,1/0.9,1/3.,0,0,1/2.9,1/0.8,1/0.8,1/0.8,1,1/0.9,1/0.9,1,1/0.9,1/0.9,1/3.4,0,0,1/2.4,1/0.8,1/0.7,1/0.8,1/0.7,1/0.8,1/0.9,1/0.9,1/0.9,1/0.9,1/2.8,0,0,1/3.2,1/1.3,1/1.1,1/1.1,1,1,1,1,1/1.1,1,1/3.6,0,0,1/2.1,1/0.7,1/0.9,1/0.8,1/0.9,1/0.9,1,1,1,1,1/2.5,0,0,1/2.3,1/0.8,1,1/0.9,1,1/0.9,1,1,1/0.9,1/1.1,1/3.4,0,0,1/2.2,1/0.9,1/1.1,1,1,1/0.9,1,1/0.9,1/0.9,1/0.8,1/3.3,0,0,1/2.,1/0.8,1/1.1,1/0.9,1/1.2,1/1.1,1/0.9,1/0.6,1/0.7,1,1/3.3,0,0,1/2.8,1/0.9,1/0.9,1/0.9,1/0.9,1/0.8,1/1.1,1/0.9,1,1,1/2.4,0,0,1/2.2,1/0.9,1/0.9,1,1,1/0.9,1/1.1,1,1/0.9,1/0.9,1/4.,0,0,1/3.3,1,1/0.9,1,1/0.8,1,1.1,1/0.8,1,1,1/3.4,0,0,1/2.7,1,1/0.8,1/0.8,1/0.8,0,1/0.8,1/0.8,1/0.7,1/0.9,1/2.7,0,0,1/2.3,1/0.9,1,1/0.9,1,1/0.8,1/0.8,1/0.9,1/0.8,1/0.9,1/2.7,0,0,1/2.3,1,1/0.8,1/0.9,1/0.9,1/0.8,1/0.8,1/0.7,1/0.9,1/1.1,1/3.1,0,0,1/3.,1/0.9,1/0.9,1/1,1/0.9,1/0.8,1/0.8,1/0.8,1/0.8,1/0.9,1/3.1,0,0,1/2.9,1/0.7,1/1.1,1/0.9,0,1/0.9,1/1.1,1/0.9,1/0.9,1,1/3.4,0,0,1/2.6,1/0.9,1/0.8,1,1/0.8,1,1/0.9,1/0.8,1/0.9,1,1/4.9,0,0,1/2.8,1/0.8,1/0.9,1,1/0.8,1,1/0.9,1/0.9,1/0.9,1,1/3.4,0,0,1/2.8,1/0.9,1,1/0.9,1/0.8,1/0.9,1,1/0.9,1,1/0.9,1/2.9,0,0,1/3.2,1/0.9,1/0.9,1/0.8,1,1/0.9,1,1/0.8,1,1/0.9,1/4.1,0,0,1/2.1,1/0.9,1/0.9,1/0.9,1/0.8,1/0.8,1/0.9,1,1,1,1/3.5,0,0,1/2.2,1/0.91,1/0.9,1,1,1/0.9,1,1/0.8,1,1,1/3.3,0,0,1/2.6,1/0.8,1/0.9,1/0.9,1/1.1,0,1/0.9,1/1.1,1/1.2,1,1/2.8,0,0,1/3.,1/0.8,1/0.9,1/0.9,1,1,1/1.1,1/1.2,1/1.1,1/1.1,1/4.9,0,0,1/1.9,1/0.8,1/0.9,1/0.8,1/0.8,1/0.8,1/0.9,1/0.7,1/0.7,1/0.9,1/2.3,0,0,1/2.2,1/0.8,1/0.8,1/0.8,1/0.9,1/0.9,1/0.9,1/0.9,1/0.9,1,1/3.6,0,0,1/2.2,1/1.6,1/1.8,1/1.3,1/2.1,1/1.7,1/1.7,1/1.8,1/1.7,1/1.2,1/2.8,0};

  for (int omnum=0; omnum<520; ++omnum){ // MW
     sncalo->setcontent(omnum,ncal[omnum]);
     // sncalo->setcontent(omnum, rate_simu[omnum]);
   }
  // sncalo->setrange(0, 3); // to force z axis range
  sncalo->draw();
  sncalo->draw1();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}



void sndisplay_calorimeter_energy_res_bdf ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");

  int om_number;
  double FWHM;
  double sigma;
  double Chi2_ndf, correction;

  double ubc_tab[520];
  memset(ubc_tab, 0, 520*sizeof(double));
  double mean_tab_simu[520];
  memset(mean_tab_simu, 0, 520*sizeof(double));
  TH1D* distrib = new TH1D("distrib", "distrib", 100, 8,25);

  TFile *file2 = new TFile("eres_1055.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("FWHM",1);
  tree2->SetBranchAddress("FWHM", &FWHM);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("Chi2_ndf",1);
  tree2->SetBranchAddress("Chi2_ndf", &Chi2_ndf);

  double means_tab[520];
  double sigmas_tab[520];
  memset (means_tab, 0, 520*sizeof(double));
  memset (sigmas_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (FWHM > 0) {
      means_tab[om_number] = FWHM;
      sigmas_tab[om_number] = sigma;
    }
  }
  int position = 0;
  TFile *file = new TFile("/eres_1055_bdf.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("FWHM",1);
  tree->SetBranchAddress("FWHM", &FWHM);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("Chi2_ndf",1);
  tree->SetBranchAddress("Chi2_ndf", &Chi2_ndf);
  tree->SetBranchStatus("position",1);
  tree->SetBranchAddress("position", &position);

  double mean_tab[520];
  double sigma_tab[520];
  memset (mean_tab, 0, 520*sizeof(double));
  memset (sigma_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (FWHM > 0 && position == 0) {
      mean_tab[om_number] = FWHM;
      sigma_tab[om_number] = sigma;
    }
  }

  double gain;
  float eres;
  int om;
  TFile *newfile = new TFile("eres_1055_bdf.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om", &om);

  std::ofstream outFile("eres_1055_bdf.txt");
  outFile << "OM" << "\t"<< "FWHM @ 1MeV" << endl;

  for (int omnum=0; omnum<520; ++omnum) // MW
    if (mean_tab[omnum] > 0){
     sncalo->setcontent(omnum, sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum]);
     if (omnum/13 != 13 && omnum/13 != 33 && omnum != 131 && omnum != 310) {
       distrib->Fill(sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum]);
     }

     eres = sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum];
     om = omnum;
     gain = mean_tab[omnum];
     outFile << om << "\t"<< eres << endl;

     Result_tree.Fill();
     // sncalo->setcontent(omnum, 2.36*sigma_tab[omnum]/mean_tab[omnum]);
   }
  // sncalo->setrange(0, 2); // to force z axis range
  distrib->Draw();
  // return;
  sncalo->draw();
  sncalo->draw1();
  outFile.close();

  newfile->cd();
  Result_tree.Write();
  newfile->Close();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}

void sndisplay_calorimeter_energy_res ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");

  int om_number;
  double mean;
  double sigma;
  double Chi2_ndf, correction;

  double ubc_tab[520];
  memset(ubc_tab, 0, 520*sizeof(double));
  double mean_tab_simu[520];
  memset(mean_tab_simu, 0, 520*sizeof(double));
  TH1D* distrib = new TH1D("distrib", "distrib", 100, 8,15);

  TFile *ubc_file = new TFile("Bi_fit/ubc_corrector.root", "READ");
  TTree* ubc_tree = (TTree*)ubc_file->Get("Result_tree");
  ubc_tree->SetBranchStatus("*",0);
  ubc_tree->SetBranchStatus("om_number",1);
  ubc_tree->SetBranchAddress("om_number", &om_number);
  ubc_tree->SetBranchStatus("correction",1);
  ubc_tree->SetBranchAddress("correction", &correction);
  ubc_tree->SetBranchStatus("mean",1);
  ubc_tree->SetBranchAddress("mean", &mean);

  for (size_t i = 0; i < 520; i++) {
    ubc_tree->GetEntry(i);
    ubc_tab[om_number] = correction;
    mean_tab_simu[om_number] = mean;
  }

  TFile *file2 = new TFile("Bi_fit/fitted_bi_gas.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("Chi2_ndf",1);
  tree2->SetBranchAddress("Chi2_ndf", &Chi2_ndf);

  double means_tab[520];
  double sigmas_tab[520];
  memset (means_tab, 0, 520*sizeof(double));
  memset (sigmas_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (mean > 0) {
      means_tab[om_number] = mean;
      sigmas_tab[om_number] = sigma;
    }
  }

  TFile *file = new TFile("Bi_fit/fitted_bi_gas_ubc.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("Chi2_ndf",1);
  tree->SetBranchAddress("Chi2_ndf", &Chi2_ndf);

  double mean_tab[520];
  double sigma_tab[520];
  memset (mean_tab, 0, 520*sizeof(double));
  memset (sigma_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mean > 0) {
      mean_tab[om_number] = mean;
      sigma_tab[om_number] = sigma;
    }
  }

  float eres;
  int om;
  TFile *newfile = new TFile("eres_simu_ubc.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om", &om);

  std::ofstream outFile("eres_simu_ubc.txt");
  outFile << "OM" << "\t"<< "FWHM @ 1MeV" << endl;
  for (int omnum=0; omnum<520; ++omnum) // MW
    if (mean_tab[omnum] > 0){
     sncalo->setcontent(omnum, sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum]);
      if (omnum/13 != 13 && omnum/13 != 33 && omnum != 131 && omnum != 310 & omnum != 74 & omnum != 80 & omnum != 98 & omnum != 119 & omnum != 136 & omnum != 146 & omnum != 461 & omnum != 369 & omnum != 318) {
        distrib->Fill(sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum]);
     }
     eres = sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum];
     om = omnum;
     outFile << om << "\t"<< eres << endl;

     Result_tree.Fill();
     // sncalo->setcontent(omnum, 2.36*sigma_tab[omnum]/mean_tab[omnum]);
  }
  distrib->Draw();

  sncalo->draw();
  sncalo->draw1();
  outFile.close();

  newfile->cd();
  Result_tree.Write();
  newfile->Close();
}

void sndisplay_calorimeter_energy_res_data ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");

  int om_number;
  double mean;
  double sigma;
  double Chi2_ndf, correction;

  double ubc_tab[520];
  memset(ubc_tab, 0, 520*sizeof(double));
  double mean_tab_simu[520];
  memset(mean_tab_simu, 0, 520*sizeof(double));
  TH1D* distrib = new TH1D("distrib", "distrib", 100, 0,20);

  TFile *file2 = new TFile("Bi_fit/fitted_bi_gas.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("mean",1);
  tree2->SetBranchAddress("mean", &mean);
  tree2->SetBranchStatus("sigma",1);
  tree2->SetBranchAddress("sigma", &sigma);
  tree2->SetBranchStatus("Chi2_ndf",1);
  tree2->SetBranchAddress("Chi2_ndf", &Chi2_ndf);

  double means_tab[520];
  double sigmas_tab[520];
  memset (means_tab, 0, 520*sizeof(double));
  memset (sigmas_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (mean > 0) {
      means_tab[om_number] = mean;
      sigmas_tab[om_number] = sigma;
    }
  }
  int position = 0;
  TFile *file = new TFile("../Bi_fit/fit_geom.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("Chi2_ndf",1);
  tree->SetBranchAddress("Chi2_ndf", &Chi2_ndf);
  tree->SetBranchStatus("position",1);
  tree->SetBranchAddress("position", &position);

  double mean_tab[520];
  double sigma_tab[520];
  memset (mean_tab, 0, 520*sizeof(double));
  memset (sigma_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (mean > 0 && position == 0) {
      // cout <<  mean_tab[omnum] << endl;
      mean_tab[om_number] = mean;
      sigma_tab[om_number] = sigma;
    }
  }

  double gain;
  float eres;
  int om;
  TFile *newfile = new TFile("eres_974_centre_ubc.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om", &om);

  std::ofstream outFile("eres_974_centre_ubc.txt");
  outFile << "OM" << "\t"<< "FWHM @ 1MeV" << endl;

  for (int omnum=0; omnum<520; ++omnum){ // MW
    if (mean_tab[omnum] > 0){
      sncalo->setcontent(omnum, sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum]);
      if (omnum/13 != 13 && omnum/13 != 33 && omnum != 131 && omnum != 310) {
        distrib->Fill(sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum]);
      }
    }
    eres = sqrt(means_tab[omnum])*235.482*sigma_tab[omnum]/mean_tab[omnum];
    if (eres > 100) cout << omnum << " -> " << eres << endl;
    om = omnum;
    gain = mean_tab[omnum];
    outFile << om << "\t"<< eres << endl;

    Result_tree.Fill();
    // sncalo->setcontent(omnum, 2.36*sigma_tab[omnum]/mean_tab[omnum]);
  }
  // sncalo->setrange(0, 2); // to force z axis range
  distrib->Draw();
  // return;
  sncalo->draw();
  sncalo->draw1();
  outFile.close();

  newfile->cd();
  Result_tree.Write();
  newfile->Close();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}

void sndisplay_calorimeter_energy_res_data_intercept ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");

  int om_number;
  double eres_1_MeV, mean;

  TFile *file = new TFile("../Bi_fit/fit_histo/fitted_bi_974_edep_AI.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("eres_1MeV",1);
  tree->SetBranchAddress("eres_1MeV", &eres_1_MeV);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);

  double mean_tab[520];
  memset (mean_tab, 0, 520*sizeof(double));
  double eres_tab[520];
  memset (eres_tab, 0, 520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (eres_1_MeV > 0) {
      // cout <<  mean_tab[omnum] << endl;
      mean_tab[om_number] = mean;
      eres_tab[om_number] = eres_1_MeV;
    }
  }

  double gain;
  float eres;
  int om;
  TFile *newfile = new TFile("eres_974_intercept.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om", &om);

  std::ofstream outFile("eres_974_intercep.txt");
  outFile << "OM" << "\t"<< "FWHM @ 1MeV" << endl;

  for (int omnum=0; omnum<520; ++omnum){ // MW
    if (eres_tab[omnum] > 0){
      sncalo->setcontent(omnum, eres_tab[omnum]);
    }
    eres = eres_tab[omnum];
    om = omnum;
    gain = mean_tab[omnum];
    outFile << om << "\t"<< eres << endl;

    Result_tree.Fill();
    // sncalo->setcontent(omnum, 2.36*sigma_tab[omnum]/mean_tab[omnum]);
  }
  // sncalo->setrange(0, 2); // to force z axis range
  // return;
  sncalo->draw();
  sncalo->draw1();
  outFile.close();

  newfile->cd();
  Result_tree.Write();
  newfile->Close();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}

void sndisplay_calorimeter_energy_res_rate ()
{
  bool with_palette = true;
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");

  int om_number;
  float eres_1_MeV, mean;

  TFile *file = new TFile("eres_974_intercept.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om",1);
  tree->SetBranchAddress("om", &om_number);
  tree->SetBranchStatus("eres",1);
  tree->SetBranchAddress("eres", &eres_1_MeV);

  double eres_tab_data[520];
  memset (eres_tab_data, 0, 520*sizeof(double));

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (eres_1_MeV > 0) {
      eres_tab_data[om_number] = eres_1_MeV;
    }
  }
  std::cout << "eres data " << eres_tab_data[10] << '\n';



  TFile *file2 = new TFile("eres_simu_ubc.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om",1);
  tree2->SetBranchAddress("om", &om_number);
  tree2->SetBranchStatus("eres",1);
  tree2->SetBranchAddress("eres", &eres_1_MeV);

  double eres_tab_simu[520];
  memset (eres_tab_simu, 0, 520*sizeof(double));

  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (eres_1_MeV > 0) {
      eres_tab_simu[om_number] = eres_1_MeV;
    }
  }
  std::cout << "eres simu " << eres_tab_simu[10] << '\n';
  for (int omnum=0; omnum<520; ++omnum){ // MW
    if (eres_tab_simu[omnum] > 0 && eres_tab_data[omnum] > 0){
      sncalo->setcontent(omnum, eres_tab_data[omnum]-eres_tab_simu[omnum]);
    }
  }
  // sncalo->setrange(0, 2); // to force z axis range
  // return;
  sncalo->draw();
  sncalo->draw1();

  // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}



// void sndisplay_calorimeter_data_comp ()
// {
//   bool with_palette = true;
//   sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
//
//   sncalo->draw_content_label("%.3f");
//
//
//   std::vector<int> *calo_type = new std::vector<int>;
//   std::vector<int> *calo_side = new std::vector<int>;
//   std::vector<int> *calo_wall = new std::vector<int>;
//   std::vector<int> *calo_column = new std::vector<int>;
//   std::vector<int> *calo_row = new std::vector<int>;
//   std::vector<int> *calo_charge = new std::vector<int>;
//   std::vector<int> *calo_ampl = new std::vector<int>;
//
//   TFile *file = new TFile("../data/snemo_run-974_udd.root", "READ");
//   TTree* tree = (TTree*)file->Get("SimData");
//   tree->SetBranchStatus("*",0);
//   tree->SetBranchStatus("digicalo.type",1);
//   tree->SetBranchAddress("digicalo.type", &calo_type);
//   tree->SetBranchStatus("digicalo.side",1);
//   tree->SetBranchAddress("digicalo.side", &calo_side);
//   tree->SetBranchStatus("digicalo.wall",1);
//   tree->SetBranchAddress("digicalo.wall", &calo_wall);
//   tree->SetBranchStatus("digicalo.column",1);
//   tree->SetBranchAddress("digicalo.column", &calo_column);
//   tree->SetBranchStatus("digicalo.row",1);
//   tree->SetBranchAddress("digicalo.row", &calo_row);
//   tree->SetBranchStatus("digicalo.charge",1);
//   tree->SetBranchAddress("digicalo.charge", &calo_charge);
//   tree->SetBranchStatus("digicalo.peakamplitude",1);
//   tree->SetBranchAddress("digicalo.peakamplitude", &calo_ampl);
//
//   double rate[712];
//   memset (rate, 0, 712*sizeof(double));
//   for (int i = 0; i < tree->GetEntries()/10; i++) {
//     tree->GetEntry(i);
//     // cout << "size = " << om_number->size() << endl;
//     for (int j = 0; j < calo_column->size(); j++) {
//       if (calo_type->at(j) == 0 && -calo_ampl->at(j) > 250) {
//         rate[calo_side->at(j)*260 + calo_column->at(j)*13 + calo_row->at(j)]++;
//         // cout << "/k" << endl;
//       }
//     }
//   }
//
//   std::vector<int> *om_number = new std::vector<int>;
//   std::vector<int> *flag_e_event = new std::vector<int>;
//
//   TFile *file2 = new TFile("../cut_974.root", "READ");
//   TTree* tree2 = (TTree*)file2->Get("Result_tree");
//   tree2->SetBranchStatus("*",0);
//   tree2->SetBranchStatus("om_number",1);
//   tree2->SetBranchAddress("om_number", &om_number);
//   tree2->SetBranchStatus("charge",1);
//   tree2->SetBranchAddress("charge", &calo_charge);
//   tree2->SetBranchStatus("amplitude",1);
//   tree2->SetBranchAddress("amplitude", &calo_ampl);
//   tree2->SetBranchStatus("flag_e_event",1);
//   tree2->SetBranchAddress("flag_e_event", &flag_e_event);
//
//   double rate_cut[712];
//   memset (rate_cut, 0, 712*sizeof(double));
//   for (int i = 0; i < tree2->GetEntries()/10; i++) {
//     tree2->GetEntry(i);
//     // cout << "size = " << om_number->size() << endl;
//     for (int j = 0; j < om_number->size(); j++) {
//       if (calo_charge->at(j) > 0 && flag_e_event->at(j) == 1 && calo_ampl->at(j) > 250) {
//         rate_cut[om_number->at(j)]++;
//         // cout << "ok" << endl;
//       }
//     }
//   }
//
//   for (size_t i = 0; i < 520; i++) {
//     std::cout << "rate cut = " << rate_cut[i] << " rate = " << rate[i] << '\n';
//   }
//
//
//
//   for (int omnum=0; omnum<520; ++omnum) // MW
//     if ( rate[omnum] > 0)
//      sncalo->setcontent(omnum, rate_cut[omnum]/rate[omnum]);
//
//   // sncalo->setrange(0, 2); // to force z axis range
//
//   sncalo->draw();
//   sncalo->draw1();
//
//   // sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
//   // sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");
//
//   // merge IT and FR canvas side by side using image magick (if installed)
//   // gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
// }


////////////////////////////////////////////////////////////////

// sndisplay_calorimeter_test_status()
// => sndisplay::calorimeter usage by filling color and/or text content

void sndisplay_calorimeter_test_status ()
{
  sndisplay::calorimeter *sncalo = new sndisplay::calorimeter ("dead_channel");

  // for (int om=0; om<712; ++om)
  //   sncalo->setcolor(om, kGreen);

  // OMs with "white photocathode PMT"
  sncalo->setmwcolor(0,  5, 9, kOrange+1); sncalo->setmwtext(0,  5, 9, "WP");
  sncalo->setmwcolor(0,  6, 2, kOrange+1); sncalo->setmwtext(0,  6, 2, "WP");
  sncalo->setmwcolor(0,  7, 7, kOrange+1); sncalo->setmwtext(0,  7, 7, "WP");
  sncalo->setmwcolor(0,  9, 2, kOrange+1); sncalo->setmwtext(0,  9, 2, "WP");
  sncalo->setmwcolor(0, 10, 6, kOrange+1); sncalo->setmwtext(0, 10, 6, "WP");
  sncalo->setmwcolor(1,  9, 5, kOrange+1); sncalo->setmwtext(1,  9, 5, "WP");
  sncalo->setmwcolor(1, 15, 6, kOrange+1); sncalo->setmwtext(1, 15, 6, "WP");

  // unfixable OMs
  sncalo->setmwcolor(0, 11, 3,   kGray); sncalo->setmwtext(0, 11, 3,   "UNFX");
  sncalo->setxwcolor(0, 1, 1, 0, kGray); sncalo->setxwtext(0, 1, 1, 0, "UNFX");
  sncalo->setxwcolor(0, 1, 0, 3, kGray); sncalo->setxwtext(0, 1, 0, 3, "UNFX");

  // new missing channels since coil installation  sncalo->setmwcolor(1, 13, 12,  kBluew); sncalo->setmwtext(1, 13, 12, "TRIP");
  sncalo->setmwcolor(0, 0, 12,   kBlue); sncalo->setmwtext(0, 0, 12,   "MISS");
  sncalo->setxwcolor(0, 1, 0, 2, kBlue); sncalo->setxwtext(0, 1, 0, 2, "MISS");

  sncalo->draw();

  sncalo->canvas_it->SaveAs("sndisplay-calorimeter-status-it.png");
  sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-status-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-status-it.png sndisplay-calorimeter-status-fr.png +append sndisplay-calorimeter-status.png");
}

////////////////////////////////////////////////////////////////

// sndisplay_calorimeter_test_omnum()
// => sndisplay::calorimeter usage to show the conversion OM_ID <=> OM_NUM

void sndisplay_calorimeter_test_omnum ()
{
  sncalo = new sndisplay::calorimeter ("omid_vs_omnum");

  for (int omnum=0; omnum<712; ++omnum)
    sncalo->settext(omnum, Form("%d", omnum));

  sncalo->draw();

  sncalo->canvas_it->SaveAs("sndisplay-calorimeter-omnum-it.png");
  sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-omnum-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  gSystem->Exec("which convert > /dev/null  && convert sndisplay-calorimeter-omnum-it.png sndisplay-calorimeter-omnum-fr.png +append sndisplay-calorimeter-omnum.png");
}

////////////////////////////////////////////////////////////////

void sndisplay_tracker_test (bool with_palette=true)
{
  sndisplay::tracker *sntracker = new sndisplay::tracker ("tracker_test", with_palette);

  // sntracker->draw_cellid_label();
  // sntracker->draw_cellnum_label();

  TRandom trand;

  for (int cellnum=0; cellnum<2034; ++cellnum)
    sntracker->setcontent(cellnum, trand.Gaus(100, 10));

  sntracker->setrange(0, 160);

  sntracker->draw();

  sntracker->canvas->SaveAs("sndisplay-tracker-test.png");
}

////////////////////////////////////////////////////////////////

void sndisplay_demonstrator_test ()
{
  sndisplay::demonstrator *sndemonstrator = new sndisplay::demonstrator ("demonstrator_test");

  const float anode_and_two_cathodes  = 1;
  const float anode_and_one_cathode   = 0.85;
  const float anode_and_no_cathode    = 0.7;
  const float two_cathodes_only       = 0.5;
  const float one_cathode_only        = 0.2;

  sndemonstrator->setomcontent(23,  1.0); // M:0.1.10 => 260*0 + 13*1 + 10 = 23
  sndemonstrator->setomcontent(276, 1.0); // M:1.1.2  => 260*1 + 13*1 + 3  = 276

  sndemonstrator->setggcontent(0,  7, 7, anode_and_two_cathodes);
  sndemonstrator->setggcontent(0,  7, 8, anode_and_one_cathode);
  sndemonstrator->setggcontent(0,  8, 6, one_cathode_only);
  sndemonstrator->setggcontent(0,  8, 7, anode_and_two_cathodes);
  sndemonstrator->setggcontent(0,  9, 2, anode_and_one_cathode);
  sndemonstrator->setggcontent(0,  9, 3, anode_and_two_cathodes);
  sndemonstrator->setggcontent(0,  9, 4, anode_and_two_cathodes);
  sndemonstrator->setggcontent(0, 10, 0, anode_and_no_cathode);
  sndemonstrator->setggcontent(0, 10, 1, anode_and_one_cathode);

  sndemonstrator->setggcontent(1,  8, 6, one_cathode_only);
  sndemonstrator->setggcontent(1,  8, 7, anode_and_one_cathode);
  sndemonstrator->setggcontent(1,  9, 3, one_cathode_only);
  sndemonstrator->setggcontent(1,  9, 4, anode_and_one_cathode);
  sndemonstrator->setggcontent(1,  9, 5, anode_and_one_cathode);
  sndemonstrator->setggcontent(1,  9, 6, anode_and_one_cathode);
  sndemonstrator->setggcontent(1,  9, 7, anode_and_one_cathode);
  sndemonstrator->setggcontent(1, 10, 0, anode_and_no_cathode);
  sndemonstrator->setggcontent(1, 10, 1, anode_and_two_cathodes);
  sndemonstrator->setggcontent(1, 10, 2, anode_and_two_cathodes);
  sndemonstrator->setggcontent(1, 10, 3, anode_and_two_cathodes);
  sndemonstrator->setggcontent(1, 11, 0, anode_and_one_cathode);

  sndemonstrator->settitle("RUN 609 // TRIGGER 9");
  sndemonstrator->draw_top();

  // disable cells in area != 0
  for (int side=0; side<2; ++side)
    for (int row=15; row<113; ++row)
      for (int layer=0; layer<9; ++layer)
	sndemonstrator->setggcolor(side, row, layer, kGray+1);

  sndemonstrator->canvas->SaveAs("sndisplay-demonstrator-test.png");
}

////////////////////////////////////////////////////////////////

// main() can be compiled with (ROOT's HistPainter library must be added!)
// g++ sndisplay.cc -o sndisplay `root-config --cflags --libs` -lHistPainter

// int main()
// {
//   sndisplay_calorimeter_test_values();
//   sndisplay_calorimeter_test_status();
//   sndisplay_calorimeter_test_omnum();

//   sndisplay_tracker_test();

//   sndisplay_demonstrator_test();

//   return 0;
// }

#endif // SNDISPLAY_CC
