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

		System.out.println("Hello");
		HipoDataSource reader = new HipoDataSource();
		H1F momentum = new H1F("momentum", "momentum", 500, 0, 10);
		H1F W_hist = new H1F("W", "W", 500, 0, 5);
		H1F Q2_hist = new H1F("Q2", "Q2", 500, 0, 10);
		H2F W_vs_Q2 = new H2F("W_vs_Q2", "W_vs_Q2", 500, 0.0, 5.0, 500, 0.0, 10.0);
		H2F mom_vs_beta = new H2F("W_vs_Q2", "W_vs_Q2", 500, 0.0, 5.0, 500, -2.5, 2.5);

		Vector3 zero_3 = new Vector3(0.0, 0.0, 0.0);
		LorentzVector p_mu = new LorentzVector();
		p_mu.setVectM(zero_3, MASS_P);
		LorentzVector e_mu = new LorentzVector(0.0, 0.0, 10.73092, 10.73092);

		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line;
			while ((line = br.readLine()) != null) {
				// process the line.

				reader.open(line);

				while (reader.hasEvent()) {
					DataEvent event = reader.getNextEvent();

					if (event.hasBank("REC::Particle")) {
						DataBank bank_rec = event.getBank("REC::Particle");

						for (int k = 0; k < bank_rec.rows(); k++) {
							int pid_rec = bank_rec.getInt("pid", k);
							float px_rec = bank_rec.getFloat("px", k);
							float py_rec = bank_rec.getFloat("py", k);
							float pz_rec = bank_rec.getFloat("pz", k);
							float beta_rec = bank_rec.getFloat("beta", k);
							float mom = (float) Math.sqrt(px_rec * px_rec + py_rec * py_rec + pz_rec * pz_rec);
							momentum.fill(mom);
							mom_vs_beta.fill(mom, beta_rec);

							if (pid_rec != 11)
								continue;
							Vector3 e_mu_3 = new Vector3(px_rec, py_rec, pz_rec);
							LorentzVector e_mu_prime = new LorentzVector();
							e_mu_prime.setVectM(e_mu_3, MASS_E);
							LorentzVector q_mu = new LorentzVector();
							q_mu.copy(e_mu);
							q_mu.sub(e_mu_prime);
							double Q2 = -q_mu.mass2();
							Q2_hist.fill(Q2);

							LorentzVector temp = new LorentzVector();
							temp.copy(p_mu);
							temp.add(q_mu);
							double W = temp.mass();
							W_hist.fill(W);
							W_vs_Q2.fill(W, Q2);

						}
					}
				}
			}
		}
		System.out.println("Hello");
		TCanvas can = new TCanvas("can", 1800, 1600);
		can.draw(W_vs_Q2, "kViridis");
		can.save("W_vs_Q2.png");

		TCanvas can4 = new TCanvas("can", 1800, 1600);
		can4.draw(mom_vs_beta, "kViridis");
		can4.save("mom_vs_beta.png");

		TCanvas can1 = new TCanvas("can", 1800, 1600);
		can1.draw(momentum);
		can1.save("mom.png");

		TCanvas can2 = new TCanvas("can", 1800, 1600);
		can2.draw(Q2_hist);
		can2.save("Q2.png");

		TCanvas can3 = new TCanvas("can", 1800, 1600);
		can3.draw(W_hist);
		can3.save("W.png");

		System.err.println("Ending this!");
		System.exit(255);

	}

}
