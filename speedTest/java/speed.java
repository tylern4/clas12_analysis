package example;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

public class example {
	static double MASS_E = 0.000511;
	static double MASS_P = 0.93827203;

	public static void main(String[] args) throws FileNotFoundException, IOException {
		String file = new String();
		if (args.length > 0) {
			file = args[0];
		} else {
			String current = new java.io.File(".").getCanonicalPath();
			file = current + "/file.lis";
		}

		HipoDataSource reader = new HipoDataSource();
		H1F momentum = new H1F("momentum", "momentum", 500, 0, 10);
		H2F mom_vs_beta = new H2F("W_vs_Q2", "W_vs_Q2", 500, 0.0, 5.0, 500, -2.5, 2.5);

		long startTime = System.nanoTime();
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			while ((line = br.readLine()) != null) {
				reader.open(line);
				while (reader.hasEvent()) {
					DataEvent event = reader.getNextEvent();
					if (event.hasBank("REC::Particle")) {
						DataBank bank_rec = event.getBank("REC::Particle");
						for (int k = 0; k < bank_rec.rows(); k++) {
							float px_rec = bank_rec.getFloat("px", k);
							float py_rec = bank_rec.getFloat("py", k);
							float pz_rec = bank_rec.getFloat("pz", k);
							float beta_rec = bank_rec.getFloat("beta", k);
							float mom = (float) Math.sqrt(px_rec * px_rec + py_rec * py_rec + pz_rec * pz_rec);
							momentum.fill(mom);
							mom_vs_beta.fill(mom, beta_rec);
						}
					}
				}
			}
		}
		long endTime = System.nanoTime();
		long duration = (endTime - startTime)/1000000.0;  //divide by 1000000 to get milliseconds.
		System.out.println(duration);

		TCanvas can0 = new TCanvas("can", 1800, 1600);
		can0.draw(mom_vs_beta, "kViridis");
		can0.save("mom_vs_beta.png");

		TCanvas can1 = new TCanvas("can", 1800, 1600);
		can1.draw(momentum);
		can1.save("mom.png");

		System.exit(255);

	}

}
