# Plan Sync

Refresh understanding of the current plan state.

## Trigger

`/plan-sync` or "sync plan", "refresh plan", "update plan state"

## Instructions

<plan-sync>

Sync your task list with the current plan file.

**Steps:**

1. **Find active plan**: Look for the most recently modified `.md` file in:
   - `~/.claude/plans/` (draft plans)
   - `project/.plans/active_plans/` (project plans)
   - `.plans/active_plans/` (project plans)

2. **Read the plan**: Use Read tool on the plan file

3. **Parse tasks**: Extract tasks from the `## Tasks` section
   - `- [ ]` = pending
   - `- [x]` = completed
   - Nested tasks under a parent inherit the parent's status context

4. **Sync task list**:
   - Use TaskList to see current tasks
   - Use TaskCreate for new tasks not in the list
   - Use TaskUpdate to mark completed tasks
   - Delete tasks that are no longer in the plan

5. **Report state**:
   ```
   Plan: {plan_name}
   File: {path}

   Progress: {completed}/{total} tasks

   Pending:
   - task 1
   - task 2

   Next action: {first pending task}
   ```

**If no plan found**: Report "No active plan found" and suggest `EnterPlanMode` if the user has a task to plan.

</plan-sync>
