/* ExternalProgram.java
 *
 * created: Fri Nov 28 2003
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998-2003  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ExternalProgramUtils.java,v 1.1 2004-06-09 09:44:34 tjc Exp $

 */
package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.Qualifier;

/**
 * Static utilities for the ExternalProgram class.
 **/
public class ExternalProgramUtils {
  public static ExternalProgramMonitor runJConProgram (final ExternalProgram external_program,
  						final FeatureVector features,
  							final StringVector sequence_file_names,
					   final Logger logger)
					   throws ExternalProgramException
					   
					    {
	try {
		   final uk.ac.sanger.jcon.JobControl jc = new uk.ac.sanger.jcon.JobControl();

		final Integer job_control_id =
		  Options.getOptions ().getIntegerProperty ("jcon_" + external_program.getName ()+
													"_program_id");

		final Integer min_jc_jobs =
		  Options.getOptions ().getIntegerProperty ("jcon_min_jobs");

		String jcon_template =
		  Options.getOptions ().getProperty ("jcon_" + external_program.getName () + "_template");

		final String jcon_batch_queue =
		  Options.getOptions ().getProperty ("jcon_batch_queue");

		   jcon_template =
			 jcon_template.replaceAll ("\\$database",  external_program.getProgramOptions ());

		   final uk.ac.sanger.jcon.manager.Administrator administrator = jc.createAdministrator ();

		   final uk.ac.sanger.jcon.job.Executable executable =
			 administrator.getExecutableById (job_control_id.intValue ());

		   final uk.ac.sanger.jcon.dao.StatusDAO status_dao = jc.createStatusDAO ();

		   final uk.ac.sanger.jcon.job.JobBatchImpl head_job =
			 new uk.ac.sanger.jcon.job.JobBatchImpl (administrator.getExecutableById (1));

		   head_job.setStatus (status_dao.readStatusById (uk.ac.sanger.jcon.job.Status.WAITING));
		   head_job.setOutputName ("/dev/null");
		   head_job.setInputName ("/dev/null");
		   head_job.setErrorName ("/dev/null");
		   head_job.setWorkDirectoryName (System.getProperty ("user.dir"));
		   head_job.setCommandTemplate ("$executable");
		   head_job.setQueue (jcon_batch_queue);

		   final  uk.ac.sanger.jcon.job.Job [] jobs = new  uk.ac.sanger.jcon.job.Job [features.size ()];

		   for (int i = 0 ; i < features.size () ; ++i) {
			 final  uk.ac.sanger.jcon.job.BatchJob job = 
			   new  uk.ac.sanger.jcon.job.JobBatchImpl (executable);

			 job.setStatus (status_dao.readStatusById ( uk.ac.sanger.jcon.job.Status.WAITING));
			 job.setWorkDirectoryName (System.getProperty ("user.dir"));

			 final String input_file_name = sequence_file_names.elementAt (i);
			 job.setInputName (input_file_name);
			 job.setOutputName (input_file_name + ".out");
			 job.setCommandTemplate (jcon_template);
			 job.setQueue (jcon_batch_queue);

			 head_job.add (job);
			 jobs[i] = job;
		   }

		   final uk.ac.sanger.jcon.job.Owner owner =
			 administrator.getOwnerByUserName (System.getProperty ("user.name"));

		   final uk.ac.sanger.jcon.job.Task task =
			   new uk.ac.sanger.jcon.job.TaskDefaultImpl (external_program.getName () + " from Artemis",
		external_program.getName () + " from Artemis",
												  owner, new uk.ac.sanger.jcon.job.Job [] { head_job });

		   final uk.ac.sanger.jcon.manager.TaskManager task_manager = jc.createTaskManager ();

		   final int task_id = task_manager.createTask (task);
		   task_manager.submitTask (task_id);

		   for (int i = 0 ; i < features.size () ; ++i) {
			 final Qualifier qualifier =
			   new Qualifier ("job",
							  "job: " + jobs[i].getId () + " " +
							  "task: " + task.getId () + " " +
							  "program: " + external_program.getName () + " " +
							  "options: " + external_program.getProgramOptions ());
			 features.elementAt (i).addQualifierValues (qualifier);
		   }
		   return new TaskMonitor (task, external_program.getName (), logger);
		 } catch (Exception e) {
		   throw new ExternalProgramException ("exception while running " +
		external_program.getName () + ": " +
											   e.getMessage ());
		 }
  }
}
